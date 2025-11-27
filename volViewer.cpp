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
buildInitialSlice( std::string name, Point lo, Point up, int axis )
{
  Point slo = lo;
  Point sup = up;
  sup[ axis ] = lo[ axis ];
  Point s = sup - slo + Point::diagonal(1);
  Point t = s + Point::diagonal(1);
  const int nb_faces    = s[ 0 ] * s[ 1 ] * s[ 2 ];
  const int nb_vertices = t[ 0 ] * t[ 1 ] * t[ 2 ];
  std::vector<RealPoint> positions( nb_vertices, RealPoint::zero );
  int i = ( axis == 0 ) ? 1 : 0; // first axis
  int j = ( axis == 2 ) ? 1 : 2; // second axis
  int idx = 0;
  RealPoint current;
  current[ axis ] = lo[ axis ] + 0.5;
  for ( int y = lo[ j ]; y <= up[ j ] + 1; y++ )
    {
      current[ j ] = y;
      for ( int x = lo[ i ]; x <= up[ i ] + 1; x++ )
        {
          current[ i ] = x;
          positions[ idx++ ] = current;
        }
    }
  idx = 0;
  std::vector<std::vector<size_t>> faces( nb_faces );
  for ( size_t y = lo[ j ]; y <= up[ j ]; y++ )
    for ( size_t x = lo[ i ]; x <= up[ i ]; x++ )
      faces[ idx++ ] = std::vector<size_t>( {
          t[ i ] * y + x,
          t[ i ] * y + x + 1,
          t[ i ] * ( y + 1 ) + x + 1,
          t[ i ] * ( y + 1 ) + x
        } );
  return polyscope::registerSurfaceMesh( name, positions, faces);
}
  
void
updateSlice( CountedPtr<SH3::GrayScaleImage> image,
             polyscope::SurfaceMesh* smesh, Point lo, Point up,
             int axis, int pos )
{
  Point slo = lo;
  Point sup = up;
  sup[ axis ] = pos;
  slo[ axis ] = pos;
  Point s = sup - slo + Point::diagonal(1);
  Point t = s + Point::diagonal(1);
  const int nb_faces    = s[ 0 ] * s[ 1 ] * s[ 2 ];
  const int nb_vertices = t[ 0 ] * t[ 1 ] * t[ 2 ];
  std::vector<RealPoint> positions( nb_vertices );
  int i = ( axis == 0 ) ? 1 : 0; // first axis
  int j = ( axis == 2 ) ? 1 : 2; // second axis
  int idx = 0;
  RealPoint current;
  current[ axis ] = pos + 0.5;
  for ( int y = lo[ j ]; y <= up[ j ] + 1; y++ )
    {
      current[ j ] = y;
      for ( int x = lo[ i ]; x <= up[ i ] + 1; x++ )
        {
          current[ i ] = x;
          positions[ idx++ ] = current;
        }
    }
  std::vector< double > colors( nb_faces );
  Domain D( slo, sup );
  idx = 0;
  for ( auto p : D )
    colors[ idx++ ] = (*image)(p);
  smesh->updateVertexPositions( positions );
  smesh->addFaceScalarQuantity( "image intensities", colors )
    ->setMapRange( { 0.0, 255.0 } )
    ->setEnabled( true );
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
    updateSlice( gray_scale_image, sliceXSurf, lo, up, 0, x );
  if ( slice_y && (ly != y) )
    updateSlice( gray_scale_image, sliceYSurf, lo, up, 1, y );
  if ( slice_z && (lz != z) )
    updateSlice( gray_scale_image, sliceZSurf, lo, up, 2, z );
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
  Point lo = K.lowerBound();
  Point up = K.upperBound();
  sliceXSurf = buildInitialSlice( "Slice X", lo, up, 0 );
  sliceYSurf = buildInitialSlice( "Slice Y", lo, up, 1 );
  sliceZSurf = buildInitialSlice( "Slice Z", lo, up, 2 );
  // extractIsosurface( threshold ); 
  // Give the hand to polyscope
  polyscope::state::userCallback = mycallback;
  polyscope::show();
  return 0;
}
