#include <iostream>
#include <vector>
#include <array>
#include <utility>

#include "CLI11.hpp"

#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>


#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"


using namespace DGtal;
using namespace Z3i;

// Using standard 3D digital space.
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;
typedef SurfaceMesh< Z3i::RealPoint, Z3i::RealVector >         SurfMesh;

// Global variables
KSpace                            K; // the space for the image
CountedPtr< SH3::BinaryImage >    binary_image;
CountedPtr< SH3::GrayScaleImage > gray_scale_image;
Parameters params;

// Global variables for GUI
int threshold = 127;
polyscope::SurfaceMesh* triSurf;
polyscope::SurfaceMesh* sliceXSurf;
polyscope::SurfaceMesh* sliceYSurf;
polyscope::SurfaceMesh* sliceZSurf;
int lx = -1;
int ly = -1;
int lz = -1;
int x  = 0;
int y  = 0;
int z  = 0;
bool slice_x = true;
bool slice_y = true;
bool slice_z = true;

polyscope::SurfaceMesh*
buildSlice( std::string name, Point lo, Point up )
{
  Domain domain( lo - Point::diagonal(1), up + Point::diagonal(1) );
  KSpace sliceK( lo - Point::diagonal(1), up + Point::diagonal(1), true );
  auto slice_image  = CountedPtr<SH3::BinaryImage>( new SH3::BinaryImage( domain ) );
  std::transform( domain.begin(), domain.end(),
                  slice_image->begin(),
                  [&] ( const Point& p ) {
                    return p.inf( lo ) == lo && p.sup( up ) == up;
                  } );
  params( "thresholdMin", 0 )( "thresholdMax", 255 );
  auto surface = SH3::makeDigitalSurface( slice_image, sliceK, params );
  auto primalSurface = SH3::makePrimalSurfaceMesh( surface );
  std::vector<std::vector<size_t>> faces;
  // Need to convert the faces
  for( size_t face= 0 ; face < primalSurface->nbFaces(); ++face )
    faces.push_back( primalSurface->incidentVertices( face ) );
  std::vector<RealPoint> positions = primalSurface->positions();
  // Color the faces
  auto surfels = SH3::getSurfelRange( surface );
  std::vector< double > colors( surfels.size() );
  for( size_t i = 0; i < surfels.size(); ++i )
    {
      auto    s   = surfels[ i ];  // current surfel
      auto    k   = sliceK.sOrthDir( s );  // its orthogonal direction
      auto    vox = sliceK.sDirectIncident( s, k ); // the incident (interior) voxel
      auto    p   = sliceK.sCoords( vox ); // the coordinates of the voxel
      colors[ i ] = (*gray_scale_image)( p ); // get value of voxel
    }
  auto sliceSurf = polyscope::registerSurfaceMesh( name, positions, faces);
  sliceSurf->addFaceScalarQuantity( "image intensities", colors )
    ->setMapRange( { 0.0, 255.0 } )
    ->setEnabled( true );
  return sliceSurf;
}

void extractDigitalSurface( int t )
{
  trace.beginBlock( "Extracting digital surface" );
  // Builds a thresholded image from a gray scale image
  Domain domain = gray_scale_image->domain();
  binary_image  = CountedPtr<SH3::BinaryImage>( new SH3::BinaryImage( domain ) );
  std::transform( domain.begin(), domain.end(),
                  binary_image->begin(),
                  [&] ( const Point& p ) { return (*gray_scale_image)(p) > t; } );
  auto surface = SH3::makeDigitalSurface( binary_image, K, params );
  trace.endBlock();
  trace.beginBlock( "Make primal surface" );
  auto primalSurface = SH3::makePrimalSurfaceMesh( surface );
  trace.endBlock();
  trace.beginBlock( "Register surface in polyscope" );
  std::vector<std::vector<size_t>> faces;
  // Need to convert the faces
  for( size_t face= 0 ; face < primalSurface->nbFaces(); ++face )
    faces.push_back( primalSurface->incidentVertices( face ) );
  std::vector<RealPoint> positions = primalSurface->positions();
  auto primalSurf = polyscope::registerSurfaceMesh( "Digital surface", positions, faces);
  trace.endBlock();
}

void extractIsosurface( int t )
{
  trace.beginBlock( "Extract isosurface" );
  auto tri  = SH3::makeTriangulatedSurface( gray_scale_image,
                                            params( "thresholdMin", t ) );
  trace.endBlock();
  trace.beginBlock( "Register surface in polyscope" );
  std::vector<std::vector<size_t>> faces;
  // Need to convert the faces
  for( size_t face= 0 ; face < tri->nbFaces(); ++face )
    faces.push_back( tri->verticesAroundFace( face ) );
  triSurf = polyscope::registerSurfaceMesh( "Isosurface", tri->positions(), faces);
  trace.endBlock();
}

// Polyscope GUI Callback
void mycallback()
{
  ImGui::SliderInt("Threshold", &threshold, 0, 255 ); //, "ratio = %.3f");
  if (ImGui::Button("Isosurface"))
    {
      extractIsosurface( threshold );
    }
  ImGui::SameLine();
  if (ImGui::Button("Digital surface"))
    {
      extractDigitalSurface( threshold );
    }
  Point lo = K.lowerBound();
  Point up = K.upperBound();
  ImGui::SliderInt("X", &x, lo[0], up[ 0 ] );
  ImGui::SliderInt("Y", &y, lo[1], up[ 1 ] );
  ImGui::SliderInt("Z", &z, lo[2], up[ 2 ] );
  ImGui::Checkbox("Slice X", &slice_x);
  ImGui::SameLine();
  ImGui::Checkbox("Slice Y", &slice_y);
  ImGui::SameLine();
  ImGui::Checkbox("Slice Z", &slice_z);
  if ( slice_x && (lx != x) )
    {
      Point slo = lo, sup = up;
      slo[0] = sup[0] = x;
      sliceXSurf = buildSlice( "Slice X", slo, sup );
    }
  if ( slice_y && (ly != y) )
    {
      Point slo = lo, sup = up;
      slo[1] = sup[1] = y;
      sliceYSurf = buildSlice( "Slice Y", slo, sup );
    }
  if ( slice_z && (lz != z) )
    {
      Point slo = lo, sup = up;
      slo[2] = sup[2] = z;
      sliceZSurf = buildSlice( "Slice Z", slo, sup );
    }
  lx = x;
  ly = y;
  lz = z;
  sliceXSurf->setEnabled( slice_x );
  sliceYSurf->setEnabled( slice_y );
  sliceZSurf->setEnabled( slice_z );
}

// main program
int main( int argc, char* argv[] )
{
  polyscope::init();
  
  CLI::App app{"Homotopic Thinning demo"};
  std::string filename;
  app.add_option("-i,--input,1", filename, "Input VOL file")->required()->check(CLI::ExistingFile);
  CLI11_PARSE(app,argc,argv);
  
  // Read voxel object and hands surfaces to polyscope
  params = SH3::defaultParameters()| SHG3::defaultParameters() | SHG3::parametersGeometryEstimation();
  params( "closed", 1)("surfaceComponents", "All");
  gray_scale_image  = SH3::makeGrayScaleImage( filename );  
  K = SH3::getKSpace( gray_scale_image );
  extractIsosurface( threshold ); 
  // Give the hand to polyscope
  polyscope::state::userCallback = mycallback;
  polyscope::show();
  return 0;
}
