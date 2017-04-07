#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "QEntropyCLP.h"

#include <fstream>
#include "itkScalarImageToHistogramGenerator.h"
#include "itkThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include <iostream>
#include <string>
#include <math.h>
#include "itkImageRegionIterator.h"

#include "entropy.h"
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

using namespace std;

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

template <class T>
int DoIt( int argc, char * argv[], T )
{
    PARSE_ARGS;
    string pathSegmented = "/temp/cine_";

    struct passwd *pw = getpwuid(getuid());
    string homedir = pw->pw_dir;
    string final = homedir + pathSegmented;

    string pathEntropy = "/temp/entropy/entropy_";
    string q_entropy = homedir + pathEntropy;
    string firstSlice;
    string lastSlice;
    string pathSlices = "/temp/slices.txt";
    string slicesFile = homedir + pathSlices;
    ifstream slices(slicesFile.c_str());
    if(slices.is_open()){
        getline(slices,firstSlice);
        getline(slices, lastSlice);
    }
    slices.close();
    for(int i = atoi(firstSlice.c_str()); i < (atoi(lastSlice.c_str())-1);i++){

        typedef    unsigned char InputPixelType;
        typedef    T     OutputPixelType;

        typedef itk::Image<InputPixelType,  2> InputImageType;
        typedef itk::Image<OutputPixelType, 2> OutputImageType;

        typedef itk::ImageFileReader<InputImageType>  ReaderType;
        typedef itk::ImageFileWriter<OutputImageType> WriterType;

        typedef itk::BinaryThresholdImageFilter<
                InputImageType, InputImageType>  FilterType;
        typedef itk::CastImageFilter<InputImageType, OutputImageType> CastType;

        typename ReaderType::Pointer reader = ReaderType::New();
        itk::PluginFilterWatcher watchReader(reader, "Read Volume",
                                             CLPProcessInformation);

        reader->SetFileName( inputVolume.c_str() );
        reader->Update();

        typedef itk::Image<unsigned short,2> ImageType2D;
        string typeTiff = ".tif";
        stringstream segment;
        segment<<final.c_str()<<(i+1)<<typeTiff;
        string filenameSegmented = segment.str();
        segment.str("");
        typename ReaderType::Pointer readerFixed = ReaderType::New();
        readerFixed->SetFileName(filenameSegmented);
        readerFixed->Update();

        typedef itk::ImageRegionIterator< InputImageType > ImageIteratorType;

        typedef itk::Image<unsigned char,2> ImageType;
        ImageType::Pointer imag = readerFixed->GetOutput();
        entropy ent;
        double threshold;
        threshold = ent.Execute(imag,qentropy);

        typename FilterType::Pointer filter = FilterType::New();

        filter->SetInput( readerFixed->GetOutput() );
        //filter->SetUseImageSpacing( useImageSpacing );
        filter->SetLowerThreshold(threshold);
        filter->SetUpperThreshold(255);
        filter->SetOutsideValue(0);
        filter->SetInsideValue(255);

        typename CastType::Pointer cast = CastType::New();
        cast->SetInput( filter->GetOutput() );

        typename WriterType::Pointer writer = WriterType::New();

        stringstream stringFile;
        stringFile<<q_entropy<<(i+1)<<typeTiff;

        string File = stringFile.str();
        stringFile.str("");

        writer->SetFileName(File);
        writer->SetInput( cast->GetOutput() );
        //writer->SetUseCompression(1);
        writer->Update();
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
