#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkPluginUtilities.h"
#include "itkImageFileReader.h"
#include "QEntropyCLP.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkThresholdImageFilter.h"
#include "itkAndImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImage.h"
#include "entropy.h"
#include <pwd.h>
#include <iostream>

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
    string pathSegmented = "/temp/segmentedFinal_";
    string pathSecondSegmented = "/temp/entropy/entropy_";

    struct passwd *pw = getpwuid(getuid());
    string homedir = pw->pw_dir;
    string final = homedir + pathSegmented;
    string secondFinal = homedir + pathSecondSegmented;

    string pathEntropy = "/temp/entropy/entropy_";
    string pathSecondEntropy = "/temp/entropy2/secondEntropy_";
    string q_entropy = homedir + pathEntropy;
    string second_q_entropy = homedir + pathSecondEntropy;
    string firstSlice;
    string lastSlice;
    string pathSlices = "/temp/slices.txt";
    string slicesFile = homedir + pathSlices;
    ifstream slices(slicesFile.c_str());
    if(slices.is_open()){
        getline(slices, firstSlice);
        getline(slices, lastSlice);
    }
    slices.close();
    if(iteration == "1Iteration" || iteration == "2Iteration"){
        for(int i = atoi(firstSlice.c_str()); i < (atoi(lastSlice.c_str()));i++){

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

            typename ReaderType::Pointer readerFixed = ReaderType::New();
            readerFixed->SetFileName(filenameSegmented);
            readerFixed->Update();

            typedef itk::Image<unsigned char,2> ImageType;
            ImageType::Pointer imag = readerFixed->GetOutput();
            entropy ent;
            double threshold;
            threshold = ent.Execute(imag,qentropy);

            typename FilterType::Pointer filter = FilterType::New();

            filter->SetInput( readerFixed->GetOutput() );
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
            writer->Update();
        }
    }
    if(iteration == "2Iteration"){
        for(int i = atoi(firstSlice.c_str()); i < (atoi(lastSlice.c_str()));i++){

            typedef    unsigned char InputPixelType;
            typedef    T     OutputPixelType;

            typedef itk::Image<unsigned char, 2>  ImageType;

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

            string typeTiff = ".tif";
            stringstream segment;
            segment<<secondFinal.c_str()<<(i+1)<<typeTiff;
            string filenameSegmented = segment.str();
            segment.str("");

            typename ReaderType::Pointer readerFixed = ReaderType::New();
            readerFixed->SetFileName(filenameSegmented);
            readerFixed->Update();

            stringstream segmentSecond;
            if(i<9)
                segmentSecond<<final.c_str()<<"00"<<(i+1)<<typeTiff;
            if(i>=9 && i<99)
                segmentSecond<<final.c_str()<<"0"<<(i+1)<<typeTiff;
            if(i>=99)
                segmentSecond<<final.c_str()<<(i+1)<<typeTiff;
            string filenameSecondSegmented = segmentSecond.str();
            segmentSecond.str("");

            typename ReaderType::Pointer readerSecond = ReaderType::New();
            readerSecond->SetFileName(filenameSecondSegmented);
            readerSecond->Update();

            typedef itk::AndImageFilter <ImageType> AndImageFilterType;
            AndImageFilterType::Pointer andFilter = AndImageFilterType::New();

            andFilter->SetInput(0, readerFixed->GetOutput());
            andFilter->SetInput(1, readerSecond->GetOutput());
            andFilter->Update();

            typedef itk::Image<unsigned char,2> ImageType;
            ImageType::Pointer imag = andFilter->GetOutput();
            entropy ent;
            double threshold;
            threshold = ent.Execute(imag,qentropy);

            typename FilterType::Pointer filter = FilterType::New();

            filter->SetInput( andFilter->GetOutput() );
            filter->SetLowerThreshold(threshold);
            filter->SetUpperThreshold(255);
            filter->SetOutsideValue(0);
            filter->SetInsideValue(255);

            typename CastType::Pointer cast = CastType::New();
            cast->SetInput( filter->GetOutput() );

            typename WriterType::Pointer writer = WriterType::New();

            stringstream stringFile;
            stringFile<<second_q_entropy<<(i+1)<<typeTiff;

            string File = stringFile.str();
            stringFile.str("");

            writer->SetFileName(File);
            writer->SetInput( cast->GetOutput() );
            writer->Update();
        }
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
