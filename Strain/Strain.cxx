#include "itkImageFileWriter.h"
#include "itkPluginUtilities.h"
#include "StrainCLP.h"
#include <pwd.h>
#include "mapping.h"

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

    for(int i = atoi(firstSlice.c_str()); i < (atoi(lastSlice.c_str()) - 4);i++){
        typedef    unsigned short InputPixelType;

        typedef itk::Image<InputPixelType,  2> InputImageType;

        typedef itk::ImageFileReader<InputImageType>  ReaderType;

        typename ReaderType::Pointer reader = ReaderType::New();
        itk::PluginFilterWatcher watchReader(reader, "Read Volume",
                                             CLPProcessInformation);

        reader->SetFileName( inputVolume.c_str() );
        reader->Update();

        typename ReaderType::Pointer readerFixed = ReaderType::New();
        typename ReaderType::Pointer readerMoving = ReaderType::New();

        typedef itk::Image<unsigned short,2> ImageType2D;
        string typeTiff = ".tif";
        stringstream segment;
        stringstream segmentMoving;

        segment<<final.c_str()<<(i+1)<<typeTiff;
        segmentMoving<<final.c_str()<<(i+4)<<typeTiff;
        string filenameSegmented = segment.str();
        string filenameSegmentedMoving = segmentMoving.str();
        segment.str("");
        segmentMoving.str("");

        readerFixed->SetFileName(filenameSegmented);
        readerFixed->Update();

        readerMoving->SetFileName(filenameSegmentedMoving);
        readerMoving->Update();

        ImageType2D::Pointer imagFixed = readerFixed->GetOutput();
        ImageType2D::Pointer imagMoving = readerMoving->GetOutput();
        Mapping map;
        map.calcMapping(imagFixed, imagMoving, i);
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
