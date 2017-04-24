#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"
#include "itkCastImageFilter.h"

#include "FractalDimensionCLP.h"
#include "fracdimension.h"
#include "itkGradientMagnitudeImageFilter.h"
#include <iostream>

#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
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
        string firstSlice;
        string lastSlice;
        string pathSegmented;
        string pathBoxCounting;

        if(edge == "Epicardium"){
            pathSegmented = "/temp/segmentedFinal_";
            pathBoxCounting = "/temp/boxCountingEpicardium2D.txt";
        }
        else{
            pathSegmented = "/temp/segmented_";
            pathBoxCounting = "/temp/boxCountingEndocardium2D.txt";
        }

        string pathSlices = "/temp/slices.txt";

        struct passwd *pw = getpwuid(getuid());
        string homedir = pw->pw_dir;
        string final = homedir + pathSegmented;
        string finalBoxCounting = homedir + pathBoxCounting;
        ofstream boxCounting2D(finalBoxCounting.c_str());
        string slicesFile = homedir + pathSlices;
        ifstream slices(slicesFile.c_str());
        if(slices.is_open()){
            getline(slices,firstSlice);
            getline(slices, lastSlice);
        }
        slices.close();
        if (boxCounting2D.is_open())
        {
            for(int i = atoi(firstSlice.c_str()); i < atoi(lastSlice.c_str());i++){
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
                    segment<<final.c_str()<<"00"<<(i+1)<<typeTiff;
                if(i>=9 && i<99)
                    segment<<final.c_str()<<"0"<<(i+1)<<typeTiff;
                if(i>=99)
                    segment<<final.c_str()<<(i+1)<<typeTiff;
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

                ImageType2D::Pointer imag = readerS->GetOutput();

                fracdimension dim;
                double dimens;
                dimens = dim.GetBoxCountingDimension2D(imag);
                //cout<<"Dimension: "<<dimens<<endl;
                boxCounting2D<<"Image: "<<i+1<<" :"<<dimens<<endl;
            }
            boxCounting2D.close();
        }
    }

    if(dimension == "Minkowski 2D"){
        string firstSlice;
        string lastSlice;
        string pathSegmented;
        string pathMinkowski;

        if(edge == "Epicardium"){
            pathSegmented = "/temp/segmentedFinal_";
            pathMinkowski = "/temp/minkowskiEpicardium2D.txt";
        }
        else{
            pathSegmented = "/temp/segmented_";
            pathMinkowski = "/temp/minkowskiEndocardium2D.txt";
        }

        struct passwd *pw = getpwuid(getuid());
        string homedir = pw->pw_dir;
        string final = homedir + pathSegmented;
        string finalMinkowski = homedir + pathMinkowski;
        ofstream minkowski2D(finalMinkowski.c_str());
        string pathSlices = "/temp/slices.txt";
        string slicesFile = homedir + pathSlices;
                ifstream slices(slicesFile.c_str());
                if(slices.is_open()){
                    getline(slices,firstSlice);
                    getline(slices, lastSlice);
                }
                slices.close();
        if (minkowski2D.is_open())
        {
            for(int i = atoi(firstSlice.c_str()); i < atoi(lastSlice.c_str());i++){
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
                    segment<<final.c_str()<<"00"<<(i+1)<<typeTiff;
                if(i>=9 && i<99)
                    segment<<final.c_str()<<"0"<<(i+1)<<typeTiff;
                if(i>=99)
                    segment<<final.c_str()<<(i+1)<<typeTiff;
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

                ImageType2D::Pointer imag = readerS->GetOutput();

                fracdimension dim;
                double dimens;
                dimens = dim.GetMinkowskiDimension2D(imag);
                //cout<<"Dimension: "<<dimens<<endl;
                minkowski2D<<"Image: "<<i+1<<" :"<<dimens<<endl;
            }
            minkowski2D.close();
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

    if(dimension == "Differencial Box Counting 2D"){
        string firstSlice;
        string lastSlice;

        string pathSegmented;
        string pathBoxCounting;

        if(edge == "Epicardium"){
            pathSegmented = "/temp/segmentedFinal_";
            pathBoxCounting = "/temp/diffboxCountingEpicardium2D.txt";
        }
        else{
            pathSegmented = "/temp/segmented_";
            pathBoxCounting = "/temp/diffboxCountingEndocardium2D.txt";
        }

        struct passwd *pw = getpwuid(getuid());
        string homedir = pw->pw_dir;
        string final = homedir + pathSegmented;
        string finalDifferencialBoxCounting = homedir + pathBoxCounting;
        ofstream boxCounting2D(finalDifferencialBoxCounting.c_str());
        string pathSlices = "/temp/slices.txt";
        string slicesFile = homedir + pathSlices;
                ifstream slices(slicesFile.c_str());
                if(slices.is_open()){
                    getline(slices,firstSlice);
                    getline(slices, lastSlice);
                }
                slices.close();
        if (boxCounting2D.is_open())
        {
            for(int i = atoi(firstSlice.c_str()); i < atoi(lastSlice.c_str());i++){
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
                string typeTiff = ".png";
                stringstream segment;
                if(i<9)
                    segment<<final.c_str()<<"00"<<(i+1)<<typeTiff;
                if(i>=9 && i<99)
                    segment<<final.c_str()<<"0"<<(i+1)<<typeTiff;
                if(i>=99)
                    segment<<final.c_str()<<(i+1)<<typeTiff;
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
                dimens = dim.GetDBCDimension2D(imag);
                //cout<<"Dimension: "<<dimens<<endl;
                boxCounting2D<<"Image: "<<i+1<<" :"<<dimens<<endl;
            }
            boxCounting2D.close();
        }
    }

    if(dimension == "Differencial Box Counting 3D"){
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
        dimens = dim.GetDBCDimension3D(imag);
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
