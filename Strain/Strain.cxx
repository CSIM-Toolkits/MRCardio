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
    string pathVector = "/temp/strain/vector_";

    string pathMyocardiumSegmented = "/temp/segmentedArea_";

    struct passwd *pw = getpwuid(getuid());
    string homedir = pw->pw_dir;
    string final = homedir + pathSegmented;
    string segmented = homedir + pathMyocardiumSegmented;
    string vector = homedir + pathVector;

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

    for(int i = (atoi(firstSlice.c_str())+1); i < (atoi(lastSlice.c_str()));i++){
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
        typename ReaderType::Pointer readerSegmented = ReaderType::New();

        typedef itk::Image<unsigned short,2> ImageType2D;
        typedef itk::VectorImage<float,2> VectorImageType;
        string typeTiff = ".tif";
        stringstream segment;
        stringstream segmentMoving;


        segment<<final.c_str()<<i<<typeTiff;
        segmentMoving<<final.c_str()<<(i+1)<<typeTiff;
        string filenameSegmented = segment.str();
        string filenameSegmentedMoving = segmentMoving.str();
        segment.str("");
        segmentMoving.str("");

        stringstream stringFileSegmented;
        if(i<=9)
            stringFileSegmented<<segmented<<"00"<<(i)<<typeTiff;
        if(i>9 && i<=99)
            stringFileSegmented<<segmented<<"0"<<(i)<<typeTiff;
        if(i>99)
            stringFileSegmented<<segmented<<(i)<<typeTiff;
        string segmentedFile = stringFileSegmented.str();
        stringFileSegmented.str("");

        readerFixed->SetFileName(filenameSegmented);
        readerFixed->Update();

        readerMoving->SetFileName(filenameSegmentedMoving);
        readerMoving->Update();

        readerSegmented->SetFileName(segmentedFile);
        readerSegmented->Update();

        ImageType2D::Pointer imagFixed = readerFixed->GetOutput();
        ImageType2D::Pointer imagMoving = readerMoving->GetOutput();
        ImageType2D::Pointer imagSegmented = readerSegmented->GetOutput();
        Mapping map;
        map.calcMapping(imagFixed, imagMoving, imagSegmented, i);

        string typeMha = ".mha";
        stringstream vectorFile;
        vectorFile<<vector.c_str()<<i<<typeMha;
        string filenameVector = vectorFile.str();
        vectorFile.str("");

        typedef itk::ImageFileReader< VectorImageType >  VectorReaderType;
        VectorReaderType::Pointer vectorReader = VectorReaderType::New();
        vectorReader->SetFileName(filenameVector);
        vectorReader->Update();

        VectorImageType::Pointer vectorImage = vectorReader->GetOutput();

        map.calcMagnitude(vectorImage, i);
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
