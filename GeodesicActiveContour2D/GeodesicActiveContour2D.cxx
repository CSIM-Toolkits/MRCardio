
#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"
#include "segment.h"

#include <iostream>

#include "GeodesicActiveContour2DCLP.h"
#include "string"
#include "qstring.h"
#include "QString"
#include "convertvolume.h"
// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
using namespace std;
namespace
{

template <class T>
int DoIt( int argc, char * argv[], T )
{
    PARSE_ARGS;
    typedef    float InputPixelType;
    typedef    unsigned char     OutputPixelType;

    typedef itk::Image<InputPixelType,  3> InputImageType;
    typedef itk::Image<OutputPixelType, 3> OutputImageType;

    typedef itk::ImageFileReader<InputImageType>  ReaderType;
    typedef itk::ImageFileWriter<OutputImageType> WriterType;


    typename ReaderType::Pointer reader = ReaderType::New();
    itk::PluginFilterWatcher watchReader(reader, "Read Volume",
                                         CLPProcessInformation);
    const char* volume = inputVolume.c_str();
    const char* volume_out = outputVolume.c_str();
    //volume = inputVolume.c_str();
    reader->SetFileName( inputVolume.c_str() );
    reader->Update();
    int first;
    int last;
    string line;
    string line2;
    ifstream myfile ("/home/gustavo/temp/slices.txt");
    if(myfile.is_open()){
        getline(myfile,line);
        first = atoi(line.c_str());
        getline(myfile,line2);
        last = atoi(line2.c_str());
    }
    myfile.close();

    //ImageType::Pointer imag_out = ImageType::New();
    typedef itk::Image<unsigned char,3> ImageType;
    ImageType::Pointer imag_out = ImageType::New();
    segment seg;
    if(axis == "Eixo Curto"){
        seg.InternalEC(first,last,sigma,sigma_min,sigma_max,propagation,curvature,advection,rms,iterations,timestep,it,conductance,alpha,beta,distance);
        seg.MyocardiumEC(first,last,sigma,sigma_min,sigma_max,propagation,curvature,advection,rms,iterations,timestep,it,conductance,alpha,beta,distance);
    }
    if(axis == "Eixo Longo Vertical"){

    }
    if(axis == "Eixo Longo Horizontal"){

    }

    convertVolume convert;
    imag_out = convert.Convert(first,last);
    typename WriterType::Pointer writer = WriterType::New();
    itk::PluginFilterWatcher watchWriter(writer, "Write Volume",
                                         CLPProcessInformation);
    writer->SetFileName( outputVolume.c_str() );
    writer->SetInput(imag_out);
    writer->SetUseCompression(1);
    writer->Update();

    return EXIT_SUCCESS;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType(inputVolume, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch( componentType )
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0) );
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<char>(0) );
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( argc, argv, static_cast<unsigned short>(0) );
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( argc, argv, static_cast<unsigned int>(0) );
        break;
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0) );
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
