/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
/**
 * @file
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systemes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 * @date 2021/09/02
 *
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 * @date 2025/12/04
 *
 * This file is part of the DGtal library.
 */
#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include <DGtal/shapes/SurfaceMesh.h>
#include <DGtal/geometry/surfaces/DigitalSurfaceRegularization.h>
#include <DGtal/dec/PolygonalCalculus.h>

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>

using namespace DGtal;
using namespace Z3i;

// Using standard 3D digital space.
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;
// The following typedefs are useful
typedef SurfaceMesh< RealPoint, RealVector >  SurfMesh;
typedef SurfMesh::Vertices                    Vertices;
typedef SurfMesh::RealPoint                   RealPoint;
typedef SurfMesh::Face                   Face;
typedef SurfMesh::Vertex                  Vertex;
typedef PolygonalCalculus<SH3::RealPoint,SH3::RealVector> PolyCalculus;

// Global variables
polyscope::SurfaceMesh* psMesh;
polyscope::SurfaceMesh* psProjMesh=nullptr;
SurfMesh                surfmesh;
PolyCalculus*           ptrCalculus = nullptr;
SHG3::RealVectors       tnormals;
SHG3::RealVectors       iinormals;
CountedPtr<SH3::DigitalSurface> surface;
CountedPtr<SH3::BinaryImage>    binary_image;
Parameters              params;
std::vector<std::vector<SH3::SurfaceMesh::Vertex>> faces;
std::vector< RealPoint > centroids;

// Other global variables
std::vector<double> phiV;
float   radiusII = 3.0;
float   scale    = 0.1;
bool useCorrectedCalculus = false;

//Restriction of an ambient scalar function to vertices
double phiVertex(const Vertex v)
{
  auto p = surfmesh.position(v);
  double x,y,z;
  std::tie(x,y,z) = { p[ 0 ], p[ 1 ], p[ 2 ] };
  return  0.5*(cos(scale*x)*sin(scale*y)
               +cos( 2.*scale*z)*sin(0.5*scale*(y-x)));
}

//Restriction of an ambient scalar function to vertices
PolyCalculus::Vector phiFace(const Face f)
{
  auto vertices = surfmesh.incidentVertices(f);
  auto nf = vertices.size();
  Eigen::VectorXd ph(nf);
  size_t cpt=0;
  for(auto v: vertices)
  {
    ph(cpt) =  phiVertex(v);
    ++cpt;
  }
  return  ph;
}

//Vertex valued function for polyscope
void initPhi()
{
  phiV.clear();
  for(auto i = 0; i < surfmesh.nbVertices(); ++i)
    phiV.push_back(phiVertex(i));
  psMesh->addVertexScalarQuantity("Phi", phiV);
}

void initQuantities()
{
  if ( ptrCalculus != nullptr ) delete ptrCalculus;
  ptrCalculus = nullptr;
  if (!useCorrectedCalculus)
    {
      ptrCalculus = new PolyCalculus(surfmesh);
      psProjMesh  = psMesh;
    }
  else
  {
    //Using the projection embedder
    ptrCalculus = new PolyCalculus(surfmesh);
    functors::EmbedderFromNormalVectors<Z3i::RealPoint, Z3i::RealVector> embedderFromNormals(iinormals,surfmesh);
    ptrCalculus->setEmbedder( embedderFromNormals );
    std::vector< PolyCalculus::Real3dVector > projPos;
    std::vector< std::vector< std::size_t > > projFaces;
    std::vector< double > projPhi;
    std::size_t idx = 0;
    for ( auto f = 0; f < faces.size(); f++ )
      {
        std::vector< std::size_t > vertices { idx, idx + 1, idx + 2, idx + 3 };
        projFaces.push_back( vertices );
        RealPoint c = RealPoint::zero;
        for ( auto v : faces[ f ] )
          {
            projPhi.push_back( phiVertex( v ) );
            auto ppos = embedderFromNormals( f, v );
            projPos.push_back( ppos );
            c += ppos;
          }
        c   /= 4;
        projPos[ idx++ ] += centroids[ f ] - c;
        projPos[ idx++ ] += centroids[ f ] - c;
        projPos[ idx++ ] += centroids[ f ] - c;
        projPos[ idx++ ] += centroids[ f ] - c;
        //idx += 4;
      }
    psProjMesh = polyscope::registerSurfaceMesh("Corrected surface", projPos, projFaces);
    psProjMesh->addVertexScalarQuantity("Phi", projPhi);
  }

  PolyCalculus& calculus = *ptrCalculus;

  std::vector<PolyCalculus::Real3dVector> gradients;
  std::vector<PolyCalculus::Vector> cogradients;
  std::vector<PolyCalculus::Real3dVector> normals;
  std::vector<PolyCalculus::Real3dVector> vectorArea;
  std::vector<PolyCalculus::Real3dPoint> centroids;
  std::vector<double> faceArea;

  for(auto f=0; f < surfmesh.nbFaces(); ++f)
  {
    PolyCalculus::Vector ph = phiFace(f);
    PolyCalculus::Vector grad = calculus.gradient(f) * ph;
    PolyCalculus::Real3dVector G( grad[0], grad[1], grad[2] );
    // Fix length by projecting onto tangent plane
    // G *= tnormals[ f ].dot( iinormals[ f ] );
    gradients.push_back( G * calculus.faceArea(f) );
    PolyCalculus::Vector cograd =  calculus.coGradient(f) * ph;
    cogradients.push_back( cograd );
    normals.push_back(calculus.faceNormalAsDGtalVector(f));
    
    auto vA = calculus.vectorArea(f);
    vectorArea.push_back({vA(0) , vA(1), vA(2)});
    
    faceArea.push_back( calculus.faceArea(f));
    
    // centroids.push_back( calculus.centroidAsDGtalPoint(f) );
  }
  
  psMesh->addFaceVectorQuantity("Gradients", gradients);
  psMesh->addFaceVectorQuantity("co-Gradients", cogradients);
  psMesh->addFaceVectorQuantity("Normals", normals);
  psMesh->addFaceScalarQuantity("Face area", faceArea);
  psMesh->addFaceVectorQuantity("Vector area", vectorArea);

  psProjMesh->addFaceVectorQuantity("Gradients", gradients);
  psProjMesh->addFaceVectorQuantity("co-Gradients", cogradients);
  psProjMesh->addFaceVectorQuantity("Normals", normals);
  psProjMesh->addFaceScalarQuantity("Face area", faceArea);
  psProjMesh->addFaceVectorQuantity("Vector area", vectorArea);
  
  //polyscope::registerPointCloud("Centroids", centroids);
}


void myCallback()
{
  ImGui::Checkbox( "Use corrected calculus", &useCorrectedCalculus );
  ImGui::SliderFloat("Phi scale", &scale, 0., 1.);
  if (ImGui::Button("Init phi"))
    initPhi();
  
  if (ImGui::Button("Compute quantities"))
    initQuantities();

  
  ImGui::SliderFloat("II radius", &radiusII , 0.,10.);
  if (ImGui::Button("Compute II normals"))
    {
      params("r-radius", (double) radiusII);
      auto surfels   = SH3::getSurfelRange( surface, params );
      iinormals = SHG3::getIINormalVectors( binary_image, surfels, params );
      trace.info()<<iinormals.size()<<std::endl;
      psMesh->addFaceVectorQuantity("II normals", iinormals);
    }
}

int main()
{
  params = SH3::defaultParameters() | SHG3::defaultParameters() |  SHG3::parametersGeometryEstimation();
  
  auto h=.3 ; //gridstep
  params( "polynomial", "goursat" )( "gridstep", h );
  auto implicit_shape  = SH3::makeImplicitShape3D  ( params );
  auto digitized_shape = SH3::makeDigitizedImplicitShape3D( implicit_shape, params );
  auto K               = SH3::getKSpace( params );
  binary_image         = SH3::makeBinaryImage( digitized_shape, params );
  surface              = SH3::makeDigitalSurface( binary_image, K, params );
  SH3::Cell2Index c2i;
  auto primalSurface   = SH3::makePrimalSurfaceMesh(c2i, surface);
  
  // Convert faces to appropriate indexed format
  for(auto face= 0 ; face < primalSurface->nbFaces(); ++face)
    faces.push_back(primalSurface->incidentVertices( face ));
  
  //Recasting to vector of vertices
  auto positions = primalSurface->positions();

  surfmesh = SurfMesh(positions.begin(),
                      positions.end(),
                      faces.begin(),
                      faces.end());
  centroids.resize( surfmesh.nbFaces() );
  for( auto f = 0; f < surfmesh.nbFaces(); ++f )
    centroids[ f ] = surfmesh.faceCentroid( f );
  
  // Initialize polyscope
  polyscope::init();
  
  psMesh = polyscope::registerSurfaceMesh("digital surface", positions, faces);

  params("r-radius", (double) radiusII);
  auto surfels   = SH3::getSurfelRange( surface, params );
  tnormals  = SHG3::getTrivialNormalVectors( K, surfels );
  iinormals = SHG3::getIINormalVectors( binary_image, surfels, params );
  trace.info()<<iinormals.size()<<std::endl;
  psMesh->addFaceVectorQuantity("II normals", iinormals);
  
  // Set the callback function
  polyscope::state::userCallback = myCallback;
  polyscope::show();
  return EXIT_SUCCESS;
}
