#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "SampleEntropyCLP.h"
#include "sampen.h"
#include "utils.h"

#include "iostream"

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
int DoIt( int argc, char * argv[], T)
{
    PARSE_ARGS;
    if(dimension == "Sample Entropy 2D"){

        string pathSegmented = "/temp/segmentedFinal_";
        string pathSampleEntropy = "/temp/sampleEntropy2D.txt";

        struct passwd *pw = getpwuid(getuid());
        string homedir = pw->pw_dir;
        string final = homedir + pathSegmented;
        string finalSampleEntropy = homedir + pathSampleEntropy;
        ofstream sampleEntropy(finalSampleEntropy.c_str());
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
        if (sampleEntropy.is_open())
        {
            for(int i = atoi(firstSlice.c_str()); i < atoi(lastSlice.c_str());i++){
                typedef    unsigned char InputPixelType;

                typedef itk::Image<InputPixelType,  2> InputImageType;

                typedef itk::ImageFileReader<InputImageType>  ReaderType;

                typename ReaderType::Pointer reader = ReaderType::New();
                itk::PluginFilterWatcher watchReader(reader, "Read Volume",
                                                     CLPProcessInformation);

                reader->SetFileName( inputVolume.c_str() );
                reader->Update();

                typename ReaderType::Pointer readerS = ReaderType::New();

                reader->SetFileName( inputVolume.c_str() );
                reader->Update();
                typedef itk::Image<unsigned char,2> ImageType2D;
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

                ImageType2D::Pointer imag = readerS->GetOutput();

                SampEn sampleEntropy2D;
                Utils stand;
                double std;
                std = stand.GetStd(imag);
                double result = sampleEntropy2D.calcSampleEn2D(imag, m , (r*std));
                //cout<<"Dimension: "<<dimens<<endl;
                sampleEntropy<<"Image: "<<i<<" :"<<result<<endl;
            }
            sampleEntropy.close();
        }
    }
    if(dimension == "Sample Entropy 3D"){

        string pathSegmented = "/temp/segmentedFinal_";
        string pathSampleEntropy3D = "/temp/sampleEntropy3D.txt";

        struct passwd *pw = getpwuid(getuid());
        string homedir = pw->pw_dir;
        string final = homedir + pathSegmented;
        string finalSampleEntropy3D = homedir + pathSampleEntropy3D;
        ofstream sampleEntropy(finalSampleEntropy3D.c_str());

        if (sampleEntropy.is_open())
        {
            typedef    unsigned char InputPixelType;

            typedef itk::Image<InputPixelType,  3> InputImageType;

            typedef itk::ImageFileReader<InputImageType>  ReaderType;

            typename ReaderType::Pointer reader = ReaderType::New();
            itk::PluginFilterWatcher watchReader(reader, "Read Volume",
                                                 CLPProcessInformation);

            reader->SetFileName( inputVolume.c_str() );
            reader->Update();

            reader->SetFileName( inputVolume.c_str() );
            reader->Update();
            typedef itk::Image<unsigned char,3> ImageType3D;

            ImageType3D::Pointer imag = reader->GetOutput();

            SampEn sampleEntropy3D;
            double result = sampleEntropy3D.calcSampleEn3D(imag, m , r);
            //cout<<"Dimension: "<<dimens<<endl;
            sampleEntropy<<"Image: "<<result<<endl;
            sampleEntropy.close();
        }
    }
    return EXIT_SUCCESS;
} // end of anonymous namespace

}

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
