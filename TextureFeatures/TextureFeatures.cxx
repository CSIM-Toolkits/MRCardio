#include "itkImageFileWriter.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkPluginUtilities.h"
#include "TextureFeaturesCLP.h"
#include "haralick.h"
#include "extractfeatures.h"
#include <pwd.h>

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

    if(dimension == "Haralick Features 2D"){
        typedef float InputPixelType;

        const unsigned int Dimension = 2;

        string pathSegmented = "/temp/segmentedFinal_";
        string pathHaralick = "/temp/haralick2D.txt";

        struct passwd *pw = getpwuid(getuid());
        string homedir = pw->pw_dir;
        string final = homedir + pathSegmented;
        string finalHaralick = homedir + pathHaralick;
        ofstream haralick2D(finalHaralick.c_str());
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
        if (haralick2D.is_open())
        {
            for(int i = atoi(firstSlice.c_str()); i < atoi(lastSlice.c_str());i++){

                typedef itk::Image<InputPixelType,  Dimension> InputImageType;

                typedef itk::ImageFileReader<InputImageType>  ReaderType;

                typename ReaderType::Pointer reader = ReaderType::New();
                typename ReaderType::Pointer readerS = ReaderType::New();

                reader->SetFileName( inputVolume.c_str() );
                reader->Update();
                typedef itk::Image<float,2> ImageType;
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
                ImageType::Pointer imag = readerS->GetOutput();
                ExtractFeatures hac;

                typedef itk::Image<float, 2> InternalImageType;
                typedef itk::Neighborhood<float, 2> NeighborhoodType;
                NeighborhoodType neighborhood;
                typedef InternalImageType::OffsetType OffsetType;
                neighborhood.SetRadius(1);
                unsigned int centerIndex = neighborhood.GetCenterNeighborhoodIndex();
                OffsetType offset;
                for ( unsigned int d = 0; d < centerIndex; d++ )
                {
                    double entropy;
                    double energy;
                    double correlation;
                    double inertia;
                    double haralickCorrelation;
                    double inverseDifferenceMoment;
                    double clusterProminence;
                    double clusterShade;

                    offset = neighborhood.GetOffset(d);
                    hac.Extract(offset,imag, &entropy, &energy, &correlation, &inertia, &haralickCorrelation,
                                &inverseDifferenceMoment, &clusterProminence, &clusterShade);
                    haralick2D<<"Image: "<<i<<endl;
                    haralick2D<<"Position: "<<d+1<<endl;
                    haralick2D<<"________________"<<endl;
                    haralick2D<<endl;
                    haralick2D<<"Entropy: "<<entropy<<endl;
                    haralick2D<<"Energy: "<<energy<<endl;
                    haralick2D<<"Correlation: "<<correlation<<endl;
                    haralick2D<<"Inertia: "<<inertia<<endl;
                    haralick2D<<"HaralickCorrelation: "<<haralickCorrelation<<endl;
                    haralick2D<<"InverseDifferenceMoment: "<<inverseDifferenceMoment<<endl;
                    haralick2D<<"ClusterProminence: "<<clusterProminence<<endl;
                    haralick2D<<"ClusterShade: "<<clusterShade<<endl;
                    haralick2D<<"------------------------------------------"<<endl;
                }

            }
            haralick2D.close();
        }
    }
    if(dimension == "Haralick Features 3D"){

        string pathHaralick = "/temp/haralick3D.txt";

        struct passwd *pw = getpwuid(getuid());
        string homedir = pw->pw_dir;
        string finalHaralick = homedir + pathHaralick;
        ofstream haralick3D(finalHaralick.c_str());
        if (haralick3D.is_open())
        {
            double entropy;
            double energy;
            double correlation;
            double inertia;
            double haralickCorrelation;
            double inverseDifferenceMoment;
            double clusterProminence;
            double clusterShade;

            typedef float InputPixelType3D;

            const unsigned int Dimension3D = 3;

            typedef itk::Image<InputPixelType3D,  Dimension3D> InputImageType3D;

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
            for ( unsigned int d = 0; d < centerIndex3D; d++ )
            {
                offset3D = neighborhood3D.GetOffset(d);
                hac3D.Extract3D(offset3D,imag3D, &entropy, &energy, &correlation, &inertia, &haralickCorrelation,
                                &inverseDifferenceMoment, &clusterProminence, &clusterShade);

                haralick3D<<"Position: "<<d+1<<endl;
                haralick3D<<"________________"<<endl;
                haralick3D<<endl;
                haralick3D<<"Entropy: "<<entropy<<endl;
                haralick3D<<"Energy: "<<energy<<endl;
                haralick3D<<"Correlation: "<<correlation<<endl;
                haralick3D<<"Inertia: "<<inertia<<endl;
                haralick3D<<"HaralickCorrelation: "<<haralickCorrelation<<endl;
                haralick3D<<"InverseDifferenceMoment: "<<inverseDifferenceMoment<<endl;
                haralick3D<<"ClusterProminence: "<<clusterProminence<<endl;
                haralick3D<<"ClusterShade: "<<clusterShade<<endl;
                haralick3D<<"------------------------------------------"<<endl;
            }

        }
        haralick3D.close();
    }
    if(dimension == "Run Length Features 3D"){

        string pathRunLength3D = "/temp/runLength3D.txt";

        struct passwd *pw = getpwuid(getuid());
        string homedir = pw->pw_dir;
        string finalRunLength = homedir + pathRunLength3D;
        ofstream runLength3D(finalRunLength.c_str());

        if (runLength3D.is_open())
        {
            double shortRunEmphasis;
            double longRunEmphasis;
            double greyLevelNonuniformity;
            double runLengthNonuniformity;
            double lowGrayLevelRunEmphasis;
            double highGreyLevelRunEmphasis;
            double shortRunLowGreyLevelEmphasis;
            double shortRunHighGreyLevelEmphasis;
            double longRunLowGreyLevelEmphasis;
            double longRunHighGreyLevelEmphasis;

            typedef float InputPixelType3D;

            const unsigned int Dimension3D = 3;

            typedef itk::Image<InputPixelType3D,  Dimension3D> InputImageType3D;

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
            OffsetType3D offset3D;

            hac3D.ExtractRunLength3D(imag3D, &shortRunEmphasis, &longRunEmphasis, &greyLevelNonuniformity,
                                     &runLengthNonuniformity, &lowGrayLevelRunEmphasis, &highGreyLevelRunEmphasis,
                                     &shortRunLowGreyLevelEmphasis, &shortRunHighGreyLevelEmphasis, &longRunLowGreyLevelEmphasis,
                                     &longRunHighGreyLevelEmphasis);
            runLength3D<<"ShortRunEmphasis: "<<shortRunEmphasis<<endl;
            runLength3D<<"LongRunEmphasis: "<<longRunEmphasis<<endl;
            runLength3D<<"GreyLevelNonuniformity: "<<greyLevelNonuniformity<<endl;
            runLength3D<<"RunLengthNonuniformity: "<<runLengthNonuniformity<<endl;
            runLength3D<<"LowGrayLevelRunEmphasis: "<<lowGrayLevelRunEmphasis<<endl;
            runLength3D<<"HighGreyLevelRunEmphasis: "<<highGreyLevelRunEmphasis<<endl;
            runLength3D<<"ShortRunLowGreyLevelEmphasis: "<<shortRunLowGreyLevelEmphasis<<endl;
            runLength3D<<"ShortRunHighGreyLevelEmphasis: "<<shortRunHighGreyLevelEmphasis<<endl;
            runLength3D<<"LongRunLowGreyLevelEmphasis: "<<longRunLowGreyLevelEmphasis<<endl;
            runLength3D<<"LongRunHighGreyLevelEmphasis: "<<longRunHighGreyLevelEmphasis<<endl;
            runLength3D<<"------------------------------------------"<<endl;

        }
        runLength3D.close();
    }
    if(dimension == "Run Length Features 2D"){
        typedef float InputPixelType;

        const unsigned int Dimension = 2;

        string pathSegmented = "/temp/segmentedFinal_";
        string pathRunLength = "/temp/runLength2D.txt";

        struct passwd *pw = getpwuid(getuid());
        string homedir = pw->pw_dir;
        string final = homedir + pathSegmented;
        string finalRunLength = homedir + pathRunLength;
        ofstream runLength2D(finalRunLength.c_str());
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

        if (runLength2D.is_open())
        {
            for(int i = atoi(firstSlice.c_str()); i < atoi(lastSlice.c_str());i++){

                typedef itk::Image<InputPixelType,  Dimension> InputImageType;

                typedef itk::ImageFileReader<InputImageType>  ReaderType;

                typename ReaderType::Pointer reader = ReaderType::New();
                typename ReaderType::Pointer readerS = ReaderType::New();

                reader->SetFileName( inputVolume.c_str() );
                reader->Update();
                typedef itk::Image<float,2> ImageType;
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
                ImageType::Pointer imag = readerS->GetOutput();
                ExtractFeatures hac;

                double shortRunEmphasis;
                double longRunEmphasis;
                double greyLevelNonuniformity;
                double runLengthNonuniformity;
                double lowGrayLevelRunEmphasis;
                double highGreyLevelRunEmphasis;
                double shortRunLowGreyLevelEmphasis;
                double shortRunHighGreyLevelEmphasis;
                double longRunLowGreyLevelEmphasis;
                double longRunHighGreyLevelEmphasis;

                typedef float InputPixelType2D;

                const unsigned int Dimension2D = 2;

                typedef itk::Image<InputPixelType2D,  Dimension2D> InputImageType2D;

                typedef itk::ImageFileReader<InputImageType2D>  ReaderType2D;

                typename ReaderType2D::Pointer reader2D = ReaderType2D::New();

                reader2D->SetFileName( inputVolume.c_str() );
                reader2D->Update();
                typedef itk::Image<float,2> ImageType2D;
                ImageType2D::Pointer imag2D = reader2D->GetOutput();
                ExtractFeatures hac2D;

                typedef itk::Image<float, 2> InternalImageType2D;
                typedef itk::Neighborhood<float, 3> NeighborhoodType2D;
                NeighborhoodType2D neighborhood2D;
                typedef InternalImageType2D::OffsetType OffsetType2D;
                neighborhood2D.SetRadius(1);
                OffsetType2D offset2D;

                hac2D.ExtractRunLength2D(imag2D, &shortRunEmphasis, &longRunEmphasis, &greyLevelNonuniformity,
                                         &runLengthNonuniformity, &lowGrayLevelRunEmphasis, &highGreyLevelRunEmphasis,
                                         &shortRunLowGreyLevelEmphasis, &shortRunHighGreyLevelEmphasis, &longRunLowGreyLevelEmphasis,
                                         &longRunHighGreyLevelEmphasis);
                runLength2D<<"Image: "<<i<<endl;
                runLength2D<<"________________"<<endl;
                runLength2D<<endl;
                runLength2D<<"ShortRunEmphasis: "<<shortRunEmphasis<<endl;
                runLength2D<<"LongRunEmphasis: "<<longRunEmphasis<<endl;
                runLength2D<<"GreyLevelNonuniformity: "<<greyLevelNonuniformity<<endl;
                runLength2D<<"RunLengthNonuniformity: "<<runLengthNonuniformity<<endl;
                runLength2D<<"LowGrayLevelRunEmphasis: "<<lowGrayLevelRunEmphasis<<endl;
                runLength2D<<"HighGreyLevelRunEmphasis: "<<highGreyLevelRunEmphasis<<endl;
                runLength2D<<"ShortRunLowGreyLevelEmphasis: "<<shortRunLowGreyLevelEmphasis<<endl;
                runLength2D<<"ShortRunHighGreyLevelEmphasis: "<<shortRunHighGreyLevelEmphasis<<endl;
                runLength2D<<"LongRunLowGreyLevelEmphasis: "<<longRunLowGreyLevelEmphasis<<endl;
                runLength2D<<"LongRunHighGreyLevelEmphasis: "<<longRunHighGreyLevelEmphasis<<endl;
                runLength2D<<"------------------------------------------"<<endl;

            }
            runLength2D.close();
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
