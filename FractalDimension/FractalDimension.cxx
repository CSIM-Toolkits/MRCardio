#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"
#include "itkCastImageFilter.h"

#include "FractalDimensionCLP.h"
#include "fracdimension.h"
#include "itkGradientMagnitudeImageFilter.h"
#include <iostream>

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

    if(dimension == "Box Counting 2D"){
        string pathSegmented = "/home/gustavo/temp/segmentedFinal_";
        ofstream boxCounting2D("/home/gustavo/temp/boxCounting2D.txt");

        if (boxCounting2D.is_open())
        {
            for(int i =1; i<100;i++){
                typedef    unsigned int InputPixelType;
                typedef    T     OutputPixelType;

                typedef itk::Image<InputPixelType,  2> InputImageType;
                typedef itk::Image<OutputPixelType, 2> OutputImageType;

                typedef itk::ImageFileReader<InputImageType>  ReaderType;
                typedef itk::ImageFileWriter<OutputImageType> WriterType;

                typedef itk::CastImageFilter<InputImageType, OutputImageType> CastType;

                typename ReaderType::Pointer reader = ReaderType::New();
                itk::PluginFilterWatcher watchReader(reader, "Read Volume",
                                                     CLPProcessInformation);

                reader->SetFileName( inputVolume.c_str() );
                reader->Update();

                typename ReaderType::Pointer readerS = ReaderType::New();

                reader->SetFileName( inputVolume.c_str() );
                reader->Update();
                typedef itk::Image<unsigned int,2> ImageType2D;
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

                // Setup types
                typedef itk::Image< unsigned int,  2 >  UnsignedCharImageType;
                typedef itk::Image< unsigned int,  2 >   FloatImageType;

                typedef itk::GradientMagnitudeImageFilter<
                        UnsignedCharImageType, FloatImageType >  filterType;


                // Create and setup a gradient filter
                filterType::Pointer gradientFilter = filterType::New();
                gradientFilter->SetInput( readerS->GetOutput() );
                gradientFilter->Update();

                ImageType2D::Pointer imag = gradientFilter->GetOutput();

                fracdimension dim;
                double dimens;
                dimens = dim.GetBoxCountingDimension2D(imag);
                //cout<<"Dimension: "<<dimens<<endl;
                boxCounting2D<<"Image: "<<i<<" :"<<dimens<<endl;
            }
            boxCounting2D.close();
        }
    }

    if(dimension == "Box Counting 3D"){
        typedef    unsigned int InputPixelType;
        typedef    T     OutputPixelType;

        typedef itk::Image<InputPixelType,  3> InputImageType;
        typedef itk::Image<OutputPixelType, 3> OutputImageType;

        typedef itk::ImageFileReader<InputImageType>  ReaderType;
        typedef itk::ImageFileWriter<OutputImageType> WriterType;

        typedef itk::CastImageFilter<InputImageType, OutputImageType> CastType;

        typename ReaderType::Pointer reader = ReaderType::New();
        itk::PluginFilterWatcher watchReader(reader, "Read Volume",
                                             CLPProcessInformation);

        reader->SetFileName( inputVolume.c_str() );
        reader->Update();

        typedef itk::ImageRegionIterator< InputImageType > ImageIteratorType;

        const char* volume_out = outputVolume.c_str();

        typedef itk::Image<unsigned int,3> ImageType;
        ImageType::Pointer imag = reader->GetOutput();
        fracdimension dim;
        double dimens;
        //dimens = dim.GetDBCDimension(imag);
        //dim.GradientMagnitude(imag, volume_out);
        dimens = dim.GetBoxCountingDimension3D(imag);
        cout<<"Dimension: "<<dimens<<endl;

    }

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
