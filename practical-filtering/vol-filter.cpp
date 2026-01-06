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
CountedPtr< SH3::GrayScaleImage > current_image;
CountedPtr< SH3::GrayScaleImage > gray_scale_image;
CountedPtr< SH3::GrayScaleImage > lung_image;
CountedPtr< SH3::GrayScaleImage > output_image;
CountedPtr< SH3::LightDigitalSurface > main_surface;
SH3::SurfelRange all_surfels;
Parameters params;

// Global variables for GUI
int threshold = 128;
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
bool refresh = true;
int last_img_choice = 1;
int      img_choice = 0;
int  closing_radius = 5;
int    minimum_size = 1000;

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

CountedPtr< SH3::GrayScaleImage >
fillSurface( const SH3::SurfelRange& surfels,
             bool inverse )
{
  Point lo = K.lowerBound();
  Point up = K.upperBound();
  Domain D( lo, up );
  // The image is filled with zero at the beginning
  CountedPtr< SH3::GrayScaleImage > output = SH3::makeGrayScaleImage( D );
  std::queue< Point > Q;
  for ( auto s : surfels )
    {
      auto k = K.sOrthDir( s );
      auto int_voxel = inverse ? K.sIndirectIncident( s, k ) : K.sDirectIncident( s, k );
      auto ext_voxel = inverse ? K.sDirectIncident( s, k ) : K.sIndirectIncident( s, k );
      auto int_p = K.sCoords( int_voxel );
      auto ext_p = K.sCoords( ext_voxel );
      Q.push( int_p );
      output->setValue( ext_p, 1 );
      output->setValue( int_p, 255 );
    }
  Point dv[ 3 ] = { Point( 1,0,0 ), Point( 0,1,0 ), Point( 0,0,1 ) };
  int n = 0;
  while ( ! Q.empty() )
    {
      if ( n++ % 1000 == 0 ) std::cout << Q.size() << "\n";
      auto p = Q.front();
      Q.pop();
      for ( int i = 0; i < 3; i++ )
        {
          const auto q = p - dv[ i ];
          if ( ( p[ i ] > lo[ i ] ) && ( (*output)( q ) == 0 ) )
            {
              Q.push( q );
              output->setValue( q, 255 );
            }
          const auto r = p + dv[ i ];
          if ( ( p[ i ] < up[ i ] ) && ( (*output)( r ) == 0 ) )
            {
              Q.push( r );
              output->setValue( r, 255 );
            }
        }
    }
  return output;
}


void extractIsosurface( CountedPtr< SH3::GrayScaleImage > image, int t )
{
  trace.beginBlock( "Extract isosurface" );
  auto tri  = SH3::makeTriangulatedSurface( image,
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

void extractDigitalSurface( CountedPtr< SH3::GrayScaleImage > image, int t )
{
  trace.beginBlock( "Extracting digital surface" );
  // Builds a thresholded image from a gray scale image
  Domain domain = image->domain();
  binary_image  = CountedPtr<SH3::BinaryImage>( new SH3::BinaryImage( domain ) );
  std::transform( domain.begin(), domain.end(),
                  binary_image->begin(),
                  [&] ( const Point& p ) { return (*image)(p) > t; } );
  auto surface  = SH3::makeDigitalSurface( binary_image, K, params );
  // Extracts a vector of connected surfaces.
  auto vec_surfs= SH3::makeLightDigitalSurfaces( binary_image, K, params );
  std::sort( vec_surfs.begin(), vec_surfs.end(),
             [] ( const auto& s1, const auto& s2) -> bool
             { return s1->size() > s2->size(); } );
  trace.info() << "#connected components = " << vec_surfs.size() << std::endl;
  all_surfels.clear();
  unsigned int n = 0;
  std::vector<int> V;
  for ( auto&& surf : vec_surfs ) {
    // Only process big enough components
    if ( surf->size() < minimum_size ) break;
    auto surfels = SH3::getSurfelRange( surf, params );
    for ( auto&& s : surfels )
      {
        all_surfels.push_back( s );
        V.push_back( n );
      }
    n += 1;
  }
  vec_surfs.resize( n );
  trace.info() << "#connected components = " << vec_surfs.size()
               << " after filtering" << std::endl;
  trace.endBlock();
  trace.beginBlock( "Register surface in polyscope" );
  std::vector<RealPoint>           positions;
  std::vector<std::vector<size_t>> faces;
  positions.reserve( V.size() + 1000 );
  faces.reserve( V.size() );
  for ( auto&& surf : vec_surfs )
    {
      std::size_t base_idx = positions.size();
      auto   primalSurface = SH3::makePrimalSurfaceMesh( surf );
      auto  surf_positions = primalSurface->positions();
      for ( auto && p : surf_positions ) positions.push_back( p );
      for( size_t face= 0 ; face < primalSurface->nbFaces(); ++face )
        {
          auto vertices = primalSurface->incidentVertices( face );
          for ( auto& v : vertices ) v += base_idx;
          faces.push_back( vertices );
        }
    }
  auto primalSurf = polyscope::registerSurfaceMesh( "Digital surface", positions, faces);
  primalSurf->addFaceScalarQuantity( "labels", V )
    ->setMapRange( {0.0, double(n)} )
    ->setColorMap( "rainbow" )
    ->setEnabled( true );
  trace.endBlock();
  main_surface = vec_surfs[0];
}

CountedPtr<SH3::GrayScaleImage>
makeDilation( CountedPtr<SH3::GrayScaleImage> image_ptr )
{
  const SH3::GrayScaleImage& image = *image_ptr;
  // Clone image
  CountedPtr<SH3::GrayScaleImage> output_ptr
    = CountedPtr<SH3::GrayScaleImage>( new SH3::GrayScaleImage( image ) );
  Domain D = image.domain();
  Point lo = D.lowerBound();
  Point up = D.upperBound();
  Point dv[ 3 ] = { Point( 1,0,0 ), Point( 0,1,0 ), Point( 0,0,1 ) };
  for ( auto p : D )
    {
      auto val = image( p );
      for ( int i = 0; i < 3; i++ )
        {
          if ( p[i] > lo[i] ) val = std::max( val, image( p - dv[i] ) );
          if ( p[i] < up[i] ) val = std::max( val, image( p + dv[i] ) );
        }
      output_ptr->setValue( p, val );
    }
  return output_ptr;
}

CountedPtr<SH3::GrayScaleImage>
makeErosion( CountedPtr<SH3::GrayScaleImage> image_ptr )
{
  const SH3::GrayScaleImage& image = *image_ptr;
  // Clone image
  CountedPtr<SH3::GrayScaleImage> output_ptr
    = CountedPtr<SH3::GrayScaleImage>( new SH3::GrayScaleImage( image ) );
  Domain D = image.domain();
  Point lo = D.lowerBound();
  Point up = D.upperBound();
  Point dv[ 3 ] = { Point( 1,0,0 ), Point( 0,1,0 ), Point( 0,0,1 ) };
  for ( auto p : D )
    {
      auto val = image( p );
      for ( int i = 0; i < 3; i++ )
        {
          if ( p[i] > lo[i] ) val = std::min( val, image( p - dv[i] ) );
          if ( p[i] < up[i] ) val = std::min( val, image( p + dv[i] ) );
        }
      output_ptr->setValue( p, val );
    }
  return output_ptr;
}

// Polyscope GUI Callback
void mycallback()
{
  ImGui::SliderInt("Closing radius", &closing_radius, 0, 15 ); //, "ratio = %.3f");
  Point lo = K.lowerBound();
  Point up = K.upperBound();
  if ( ImGui::Button( "Dilation" ) )
    {
      *current_image = *makeDilation( current_image );
      refresh = true;
    }
  ImGui::SameLine();
  if ( ImGui::Button( "Erosion" ) )
    {
      *current_image = *makeErosion( current_image );
      refresh = true;
    }
  ImGui::SameLine();
  if ( ImGui::Button( "Save" ) )
    {
      std::cout << "Saving <output.vol> image... ";
      bool ok = SH3::saveGrayScaleImage( current_image, "output.vol" );
      std::cout << ( ok ? "ok." : "error." ) << std::endl;
    }

  ImGui::SliderInt("Threshold", &threshold, 0, 255 ); //, "ratio = %.3f");
  ImGui::SliderInt("Minimum size", &minimum_size, 0, 10000 );
  if (ImGui::Button("Isosurface"))
    {
      extractIsosurface( current_image, threshold );
    }
  ImGui::SameLine();
  if (ImGui::Button("Digital surface"))
    {
      extractDigitalSurface( current_image, threshold );
    }
  ImGui::SameLine();
  if (ImGui::Button("Fill main surf."))
    {
      lung_image = fillSurface( SH3::getSurfelRange( main_surface, params ),
                                true );
      refresh = true;
    }
  ImGui::SameLine();
  if (ImGui::Button("Select lungs"))
    {
      for ( int i = 0; i < closing_radius; i++ )
        lung_image = makeDilation( lung_image );
      for ( int i = 0; i < closing_radius; i++ )
        lung_image = makeErosion( lung_image );
      // filter input image
      auto itOut = output_image->begin();
      for ( auto it = lung_image->begin(), ite = lung_image->end();
            it != ite; ++it )
        *itOut++ = (*it) > 128 ? *itOut : 0;
      refresh = true;
    }
  ImGui::SameLine();
  if (ImGui::Button("Fill all surf."))
    {
      output_image = fillSurface( all_surfels, false );
      refresh = true;
    }
  
  ImGui::RadioButton("Input image",  &img_choice, 0); ImGui::SameLine();
  ImGui::RadioButton("Lung image",   &img_choice, 1); ImGui::SameLine();
  ImGui::RadioButton("Output image", &img_choice, 2);
  current_image =
    img_choice == 0   ? gray_scale_image
    : img_choice == 1 ? lung_image
    :                   output_image;
  if ( last_img_choice != img_choice ) refresh = true;
  last_img_choice = img_choice;
  
  ImGui::SliderInt("X", &x, lo[0], up[ 0 ] );
  ImGui::SliderInt("Y", &y, lo[1], up[ 1 ] );
  ImGui::SliderInt("Z", &z, lo[2], up[ 2 ] );
  ImGui::Checkbox("Slice X", &slice_x);
  ImGui::SameLine();
  ImGui::Checkbox("Slice Y", &slice_y);
  ImGui::SameLine();
  ImGui::Checkbox("Slice Z", &slice_z);
  if ( refresh || ( slice_x && (lx != x) ) )
    updateSlice( current_image, sliceXSurf, lo, up, 0, x );
  if ( refresh || ( slice_y && (ly != y) ) )
    updateSlice( current_image, sliceYSurf, lo, up, 1, y );
  if ( refresh || ( slice_z && (lz != z) ) )
    updateSlice( current_image, sliceZSurf, lo, up, 2, z );
  lx = x;
  ly = y;
  lz = z;
  sliceXSurf->setEnabled( slice_x );
  sliceYSurf->setEnabled( slice_y );
  sliceZSurf->setEnabled( slice_z );
  refresh = false;
}

// main program
int main( int argc, char* argv[] )
{
  polyscope::init();
  
  CLI::App app{"Filtering practical"};
  std::string filename;
  app.add_option("-i,--input,1", filename, "Input VOL file")->required()->check(CLI::ExistingFile);
  CLI11_PARSE(app,argc,argv);
  
  // Read voxel object and hands everything to polyscope
  params = SH3::defaultParameters()| SHG3::defaultParameters() | SHG3::parametersGeometryEstimation();
  params( "closed", 1)("surfaceComponents", "All")("surfelAdjacency", 1);
  gray_scale_image  = SH3::makeGrayScaleImage( filename );
  lung_image        = gray_scale_image;
  output_image      = gray_scale_image;
  K = SH3::getKSpace( gray_scale_image );
  
  Point lo = K.lowerBound();
  Point up = K.upperBound();
  sliceXSurf = buildInitialSlice( "Slice X", lo, up, 0 );
  sliceYSurf = buildInitialSlice( "Slice Y", lo, up, 1 );
  sliceZSurf = buildInitialSlice( "Slice Z", lo, up, 2 );
  // Give the hand to polyscope
  polyscope::state::userCallback = mycallback;
  polyscope::show();
  return 0;
}
