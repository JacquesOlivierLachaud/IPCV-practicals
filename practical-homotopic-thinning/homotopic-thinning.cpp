#include <iostream>
#include <vector>
#include <array>
#include <utility>

#include "CLI11.hpp"

#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include <DGtal/shapes/SurfaceMesh.h>

#include <DGtal/topology/CubicalComplex.h>
#include <DGtal/topology/CubicalComplexFunctions.h>

#include <DGtal/topology/VoxelComplex.h>
#include <DGtal/topology/VoxelComplexFunctions.h>
#include <DGtal/topology/NeighborhoodConfigurations.h>
#include <DGtal/topology/tables/NeighborhoodTables.h>
// Distance Map
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/VoronoiMap.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"


#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"


using namespace DGtal;
using namespace Z3i;

// Using standard 3D digital space.
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;
typedef SurfaceMesh< Z3i::RealPoint, Z3i::RealVector >         SurfMesh;

KSpace K;
CountedPtr< SH3::BinaryImage > original_image;
CountedPtr< SH3::BinaryImage > binary_image;
CountedPtr< Z3i::Object26_6 >  the_object;
int persistence = 1;


CountedPtr< Z3i::Object26_6 >
makeObjectFromImage( CountedPtr< SH3::BinaryImage > image )
{
  const auto K = SH3::getKSpace( image );
  Domain domain( K.lowerBound(), K.upperBound() );
  Z3i::DigitalSet voxel_set( domain );
  for ( const auto & p : domain )
    if ( (*binary_image)( p ) ) voxel_set.insert( p );
  CountedPtr< Z3i::Object26_6 > object;
  object = CountedPtr< Z3i::Object26_6 >( new Z3i::Object26_6( dt26_6, voxel_set ) );
  object->setTable(functions::loadTable<3>(simplicity::tableSimple26_6));
  return object;
}


/// Register to polyscope the boundary surfels of a given binary image
/// \a bimage.
void registerDigitalSurface( CountedPtr< SH3::BinaryImage > bimage,
                             std::string name )
{
  auto params = SH3::defaultParameters() | SHG3::defaultParameters() |  SHG3::parametersGeometryEstimation();
  auto h=1.; //gridstep
  params( "closed", 1)("surfaceComponents", "All");
  auto K             = SH3::getKSpace( bimage );
  trace.beginBlock( "Extracting digital surface" );
  auto surface       = SH3::makeDigitalSurface( bimage, K, params );
  trace.endBlock();
  trace.beginBlock( "Make primal surface" );
  auto primalSurface = SH3::makePrimalSurfaceMesh(surface);
  trace.endBlock();
  trace.beginBlock( "Register surface in polyscope" );
  std::vector<std::vector<size_t>> faces;
  std::vector<RealPoint> positions;
  //Need to convert the faces
  for(size_t face= 0 ; face < primalSurface->nbFaces(); ++face)
    faces.push_back(primalSurface->incidentVertices( face ));
  //Recasting to vector of vertices
  positions = primalSurface->positions();
  // auto surfmesh = SurfMesh(positions.begin(),
  //                          positions.end(),
  //                          faces.begin(),
  //                          faces.end());
  auto primalSurf = polyscope::registerSurfaceMesh( name, positions, faces);
  trace.endBlock();
}

// Removes a peel of simple points onto voxel object.
bool oneStep( CountedPtr< Z3i::Object26_6 > object )
{
  DigitalSet & S = object->pointSet();
  std::queue< Point > Q;
  for ( auto&& p : S )
    if ( object->isSimple( p ) )
      Q.push( p );
  int nb_simple = 0;
  while ( ! Q.empty() )
    {
      const auto p = Q.front();
      Q.pop();
      if ( object->isSimple( p ) )
        {
          S.erase( p );
          binary_image->setValue( p, false );
          ++nb_simple;
        }
    }
  
  //Visualization
  trace.info() << "Removed " << nb_simple << " / " << S.size()
               << " points." << std::endl;
  return nb_simple == 0;
}

void criticalKernels( bool OneIsthmus = true, int persistence = 1)
{
  const bool verbose = false;
  DigitalSet & image_set = the_object->pointSet();
  // Create a VoxelComplex from the set

  using DigitalTopology = DT26_6;
  using Complex = DGtal::VoxelComplex<KSpace>;

  KSpace ks;
  Point d1( KSpace::Point::diagonal( 1 ) );
  ks.init(binary_image->domain().lowerBound() - d1 , binary_image->domain().upperBound() + d1 , true);

  trace.beginBlock("construct with table");
  Complex vc(ks);
  vc.construct(image_set, functions::loadTable(simplicity::tableSimple26_6 ));
  trace.endBlock();
  trace.beginBlock("load isthmus table");
  boost::dynamic_bitset<> isthmus_table;
  if ( ! OneIsthmus )
    isthmus_table = *functions::loadTable(isthmusicity::tableIsthmus);
  else 
    isthmus_table = *functions::loadTable(isthmusicity::tableOneIsthmus);
  auto pointMap = *functions::mapZeroPointNeighborhoodToConfigurationMask<Point>();
  trace.endBlock();

  using namespace DGtal::functions ;
  // SKEL FUNCTION:
  std::function< bool(const Complex&, const Cell&) > Skel ;
  //if (sk == "ulti") Skel = skelUltimate<Complex>;
  //else if (sk == "end") Skel = skelEnd<Complex>;
  Skel = [&isthmus_table, &pointMap](const Complex & fc,
                                     const Complex::Cell & c)
  {
    return skelWithTable(isthmus_table, pointMap, fc, c);
  };
  auto start = std::chrono::system_clock::now();

  // SELECT FUNCTION
  /*
   * Calculate distance map even if not requested:
   */
  trace.beginBlock("Create Distance Map");
  using L3Metric = ExactPredicateLpSeparableMetric<Z3i::Space, 3>;
  using DT       = DistanceTransformation<Z3i::Space, DigitalSet, L3Metric>;
  L3Metric l3;
  DT dt(binary_image->domain(), image_set, l3);
  trace.endBlock();

  std::function< std::pair<typename Complex::Cell, typename Complex::Data>(const Complex::Clique&) > Select ;
  Select =
    [&dt](const Complex::Clique & clique){
      return selectMaxValue<DT, Complex>(dt,clique);
    };
  
  trace.beginBlock("Thinning");
  Complex vc_new(ks);
  if (persistence == 0)
    vc_new = asymetricThinningScheme< Complex >( vc, Select,  Skel, verbose);
  else
    vc_new = persistenceAsymetricThinningScheme< Complex >( vc, Select,  Skel, persistence,
                                                            verbose);
  trace.endBlock();

  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::seconds> (end - start) ;
  std::cout <<"Time elapsed: " << elapsed.count() << std::endl;

  DigitalSet thin_set(binary_image->domain());
  vc_new.dumpVoxels(thin_set);
  for ( auto p : image_set )
    if ( thin_set.find( p ) == thin_set.end() )
      binary_image->setValue( p, false );
  
  registerDigitalSurface( binary_image, "Critical voxels" );
}

// Polyscope GUI Callback
void mycallback()
{
  if (ImGui::Button("Reset"))
    {
      binary_image = CountedPtr< SH3::BinaryImage >( new SH3::BinaryImage( *original_image ) );
      the_object = makeObjectFromImage( binary_image );
      registerDigitalSurface( binary_image, "Primal surface" );  
    }
  if (ImGui::Button("Run"))
    {
      oneStep( the_object );
      registerDigitalSurface( binary_image, "Thinned object" );
    }
  ImGui::SameLine();
  if (ImGui::Button("Run 10"))
    {
      for (int i = 0; i < 10; i++ )
        oneStep( the_object );
      registerDigitalSurface( binary_image, "Thinned object" );
    }
  ImGui::SameLine();
  if (ImGui::Button("Run 100"))
    {
      for (int i = 0; i < 100; i++ )
        oneStep( the_object );
      registerDigitalSurface( binary_image, "Thinned object" );
    }
  if (ImGui::Button("All screenshots"))
    {
      bool finished = false;
      while ( ! finished )
        {
          finished = oneStep( the_object );
          registerDigitalSurface( binary_image, "Thinned object" );
          polyscope::screenshot();
          polyscope::refresh();
        }
    }
  ImGui::SliderInt( "persistence", &persistence, 1, 10 );
  if (ImGui::Button("Isthmus"))
    criticalKernels( false, persistence );
  ImGui::SameLine();
  if (ImGui::Button("1-Isthmus"))
    criticalKernels( true, persistence );
}

// main program
int main( int argc, char* argv[] )
{
  polyscope::init();

  CLI::App    app{"Homotopic Thinning demo"};
  std::string filename;
  int         threshold_min = 0;
  int         threshold_max = 255;
  app.add_option("-i,--input,1", filename, "Input VOL file")->required()->check(CLI::ExistingFile);
  app.add_option("-m,--threshold_min,1", threshold_min, "Threshold min value (excluded)");
  app.add_option("-M,--threshold_max,1", threshold_max, "Threshold max value (included)");
  CLI11_PARSE(app,argc,argv);

  // Read voxel object and hands surfaces to polyscope
  auto params = SH3::defaultParameters()| SHG3::defaultParameters() | SHG3::parametersGeometryEstimation();
  params( "thresholdMin", threshold_min )( "thresholdMax", threshold_max );
  original_image = SH3::makeBinaryImage(filename, params );
  binary_image = CountedPtr< SH3::BinaryImage >( new SH3::BinaryImage( *original_image ) );
  //Visualization
  registerDigitalSurface( binary_image, "Primal surface" );  
  // registerDigitalSurface( binary_image, "Thinned object" );

  the_object = makeObjectFromImage( binary_image );
  
  // // Build object with digital topology
  // const auto K = SH3::getKSpace( binary_image );
  // Domain domain( K.lowerBound(), K.upperBound() );
  // Z3i::DigitalSet voxel_set( domain );
  // for ( const auto & p : domain )
  //   if ( (*binary_image)( p ) ) voxel_set.insert( p );
  
  // the_object = CountedPtr< Z3i::Object26_6 >( new Z3i::Object26_6( dt26_6, voxel_set ) );
  // the_object->setTable(functions::loadTable<3>(simplicity::tableSimple26_6));

  // Give the hand to polyscope
  polyscope::state::userCallback = mycallback;
  polyscope::show();
  return 0;
}
