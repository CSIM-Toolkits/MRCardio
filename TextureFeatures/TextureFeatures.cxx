#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkLabelImageToLabelMapFilter.h"

#include "itkPluginUtilities.h"

#include "TextureFeaturesCLP.h"
#include "haralick.h"
#include "iostream"
#include "extractfeatures.h"
using namespace std;

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

template <typename TPixel>
int DoIt( int argc, char * argv[], TPixel )
{
    PARSE_ARGS;

    if(dimension == "2D"){
        typedef float InputPixelType;
        typedef float OutputPixelType;

        const unsigned int Dimension = 2;
        string pathSegmented = "/home/gustavo/temp/segmentedFinal_";
        for(int i =1; i<100;i++){
            typedef itk::Image<InputPixelType,  Dimension> InputImageType;
            typedef itk::Image<OutputPixelType, Dimension> OutputImageType;

            typedef itk::ImageFileReader<InputImageType>  ReaderType;

            typename ReaderType::Pointer reader = ReaderType::New();
            typename ReaderType::Pointer readerS = ReaderType::New();

            reader->SetFileName( inputVolume.c_str() );
            reader->Update();
            typedef itk::Image<float,2> ImageType;
            string typeTiff = ".tif";
            stringstream segment;
            if(i<9)
                segment<<pathSegmented<<"00"<<(i+1)<<typeTiff;
            if(i>=9 && i<99)
                segment<<pathSegmented<<"0"<<(i+1)<<typeTiff;
            if(i>=99)
                segment<<pathSegmented<<(i+1)<<typeTiff;
            string filenameSegmented = segment.str();
            segment.str("");

            readerS->SetFileName(filenameSegmented);
            readerS->Update();
            ImageType::Pointer imag = readerS->GetOutput();
            ExtractFeatures hac;

            typedef itk::Image<float, 2> InternalImageType;
            typedef itk::Neighborhood<float, 2> NeighborhoodType;
            NeighborhoodType neighborhood;
            typedef InternalImageType::OffsetType OffsetType;
            neighborhood.SetRadius(1);
            unsigned int centerIndex = neighborhood.GetCenterNeighborhoodIndex();
            OffsetType offset;
            const char* volume = inputVolume.c_str();
            for ( unsigned int d = 0; d < centerIndex; d++ )
            {
                offset = neighborhood.GetOffset(d);
                cout<<endl;
                hac.Extract(offset,imag);
            }
        }
    }
    else{
        typedef float InputPixelType3D;
        typedef float OutputPixelType3D;

        const unsigned int Dimension3D = 3;

        typedef itk::Image<InputPixelType3D,  Dimension3D> InputImageType3D;
        typedef itk::Image<OutputPixelType3D, Dimension3D> OutputImageType3D;

        typedef itk::ImageFileReader<InputImageType3D>  ReaderType3D;

        typename ReaderType3D::Pointer reader3D = ReaderType3D::New();

        reader3D->SetFileName( inputVolume.c_str() );
        reader3D->Update();
        typedef itk::Image<float,3> ImageType3D;
        ImageType3D::Pointer imag3D = reader3D->GetOutput();
        ExtractFeatures hac3D;

        typedef itk::Image<float, 3> InternalImageType3D;
        typedef itk::Neighborhood<float, 3> NeighborhoodType3D;
        NeighborhoodType3D neighborhood3D;
        typedef InternalImageType3D::OffsetType OffsetType3D;
        neighborhood3D.SetRadius(1);
        unsigned int centerIndex3D = neighborhood3D.GetCenterNeighborhoodIndex();
        OffsetType3D offset3D;
        const char* volume3D = inputVolume.c_str();
        for ( unsigned int d = 0; d < centerIndex3D; d++ )
        {
            offset3D = neighborhood3D.GetOffset(d);
            hac3D.Extract3D(offset3D,imag3D);
        }
    }
    //hac.Extract(offset,imag);

    /*typedef itk::SmoothingRecursiveGaussianImageFilter<
    InputImageType, OutputImageType>  FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetSigma( sigma );

  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( filter->GetOutput() );
  writer->SetUseCompression(1);
  writer->Update();
  */

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
            return DoIt( argc, argv, static_cast<signed char>(0) );
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
            std::cerr << "Unknown input image pixel component type: ";
            std::cerr << itk::ImageIOBase::GetComponentTypeAsString( componentType );
            std::cerr << std::endl;
            return EXIT_FAILURE;
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
