#include "itkHoughTransformRadialVotingImageFilter.h"
// Software Guide : EndCodeSnippet

#include "itkImage.h"
#include "itkImageSource.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include <itkGradientMagnitudeImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <list>
#include "itkCastImageFilter.h"
#include "vnl/vnl_math.h"

#include "time.h"

#include "extract.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include <list>
#include "itkHoughTransform2DCirclesImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkMath.h"
#include <iostream>
#include "QString"
#include "itkImage.h"
#include "itkExtractImageFilter.h"
#include <string.h>
#include <sstream>
#include <itkCastImageFilter.h>
#include "itkRescaleIntensityImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include <fstream>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
using namespace std;

Extract::Extract()
{
    this->pw = getpwuid(getuid());

    this->homedir = this->pw->pw_dir;

    this->endocardium = "/temp/endocardium.txt";
    this->radius = "/temp/radius.txt";
    this->slices = "/temp/slices.txt";
    this->hough = "/temp/hough_";
    this->cine = "/temp/cine_";
    this->pathEndocardium = this->homedir + this->endocardium;
    this->pathRadius = this->homedir + this->radius;
    this->pathSlices = this->homedir + this->slices;
    this->pathHough = this->homedir + this->hough;
    this->pathCine = this->homedir + this->cine;

}

void Extract::Execute(int first,int last,const char* volume, int minimum, int maximum){

    typedef   unsigned char   PixelType;
    const     unsigned int    Dimension = 2;

    typedef itk::Image< PixelType, Dimension >  ImageType;

    typedef   float   floatPixelType;
    const     unsigned int    Dimension3 = 3;

    typedef itk::Image< floatPixelType, Dimension3 >  fImageType;

    typedef itk::Image<PixelType, 3> ImageType3D;
    typedef itk::Image<PixelType, 2> ImageType2D;

    //typedef    unsigned char InputPixelType;

    //typedef itk::Image<InputPixelType,  3> InputImageType;

    typedef   unsigned char   PixelType;
    typedef   float           AccumulatorPixelType;
    typedef itk::Image< PixelType, Dimension >  ImageType;
    ImageType::IndexType localIndex;
    typedef itk::Image< AccumulatorPixelType, Dimension > AccumulatorImageType;

    typedef  itk::ImageFileReader< fImageType > ReaderType3D;
    ReaderType3D::Pointer reader3D = ReaderType3D::New();
    reader3D->SetFileName(volume);
    reader3D->Update();

    typedef itk::ExtractImageFilter<ImageType3D, ImageType2D>  ExtractFilterType;
    typename ExtractFilterType::Pointer extractImg = ExtractFilterType::New();

    typedef itk::HoughTransform2DCirclesImageFilter<PixelType,
            AccumulatorPixelType> HoughTransformFilterType;
    HoughTransformFilterType::Pointer houghFilter
            = HoughTransformFilterType::New();

    int sizeX = reader3D->GetOutput()->GetRequestedRegion().GetSize()[0];
    int sizeY = reader3D->GetOutput()->GetRequestedRegion().GetSize()[1];
    //int sizeZ = reader3D->GetOutput()->GetRequestedRegion().GetSize()[2];

    typedef itk::RescaleIntensityImageFilter< fImageType, ImageType3D > RescaleFilterType;
    RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput(reader3D->GetOutput());
    rescaleFilter->SetOutputMinimum(0);
    rescaleFilter->SetOutputMaximum(255);
    rescaleFilter->Update();

    /*typedef itk::CastImageFilter< fImageType, ImageType3D > CastFilterType;
     CastFilterType::Pointer castFilter = CastFilterType::New();

     castFilter->SetInput(rescaleFilter->GetOutput());
     castFilter->Update();*/

    ofstream endocardiumFile (pathEndocardium.c_str());
    ofstream radiusFile (pathRadius.c_str());
    ofstream sliceFile (pathSlices.c_str());
    sliceFile<<first<<"\n";
    sliceFile<<last<<"\n";
    for(int i=first; i<last; i++){
        extractImg->SetInput(rescaleFilter->GetOutput());
        ImageType3D::RegionType region2;
        region2.SetSize(0, sizeX);
        region2.SetSize(1, sizeY);
        region2.SetSize(2, 0);

        region2.SetIndex(0, 0);
        //region.SetIndex(1, slice);
        region2.SetIndex(1, 0);
        region2.SetIndex(2, i);
        extractImg->SetExtractionRegion(region2);
#if ITK_VERSION_MAJOR >= 4
        extractImg->SetDirectionCollapseToIdentity(); // This is required.
#endif
        extractImg->Update();

        typedef   itk::GradientMagnitudeRecursiveGaussianImageFilter<
                ImageType,
                ImageType >  GradientFilterType;
        GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
        gradientMagnitude->SetInput(extractImg->GetOutput());
        gradientMagnitude->SetSigma(0.5);
        gradientMagnitude->Update();
        typedef   itk::SigmoidImageFilter<
                ImageType,
                ImageType >  SigmoidFilterType;
        SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();

        sigmoid->SetOutputMinimum(0.0);
        sigmoid->SetOutputMaximum(255.0);
        sigmoid->SetAlpha(-0.5);
        sigmoid->SetBeta(9.0);
        sigmoid->SetInput(gradientMagnitude->GetOutput());
        sigmoid->Update();

        ImageType::Pointer localImage = extractImg->GetOutput();
        houghFilter->SetInput( extractImg->GetOutput() );
        houghFilter->SetNumberOfCircles(1);
        houghFilter->SetMinimumRadius(minimum);
        houghFilter->SetMaximumRadius(maximum);
        //houghFilter->SetSweepAngle(5);

        houghFilter->SetSigmaGradient(1);

        //houghFilter->SetVariance(5);

        // houghFilter->SetDiscRadiusRatio( atof(argv[9]) );

        houghFilter->Update();

        AccumulatorImageType::Pointer localAccumulator = houghFilter->GetOutput();

        HoughTransformFilterType::CirclesListType circles;
        circles = houghFilter->GetCircles(2);

        typedef  unsigned char                            OutputPixelType;
        typedef  itk::Image< OutputPixelType, Dimension > OutputImageType;
        OutputImageType::Pointer  localOutputImage = OutputImageType::New();
        OutputImageType::RegionType region;
        region.SetSize(localImage->GetLargestPossibleRegion().GetSize());
        region.SetIndex(localImage->GetLargestPossibleRegion().GetIndex());
        localOutputImage->SetRegions( region );
        localOutputImage->SetOrigin(localImage->GetOrigin());
        localOutputImage->SetSpacing(localImage->GetSpacing());
        localOutputImage->Allocate(); // initializes buffer to zero
        //localOutputImage->FillBuffer(0);
        typedef HoughTransformFilterType::CirclesListType CirclesListType;
        CirclesListType::const_iterator itCircles = circles.begin();
        while( itCircles != circles.end() )
        {
            if (endocardiumFile.is_open() && radiusFile.is_open())
            {
                endocardiumFile<<(long int)(*itCircles)->GetObjectToParentTransform()->GetOffset()[0]<<' ';
                endocardiumFile<<(long int)(*itCircles)->GetObjectToParentTransform()->GetOffset()[1]<<"\n";
                radiusFile<<(long int)(*itCircles)->GetRadius()[0]<<"\n";
            }
            for(double angle = 0;angle <= 2*vnl_math::pi; angle += vnl_math::pi/60.0 )
            {
                localIndex[0] =
                        (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[0]
                        + (*itCircles)->GetRadius()[0]*std::cos(angle));
                localIndex[1] =
                        (long int)((*itCircles)->GetObjectToParentTransform()->GetOffset()[1]
                        + (*itCircles)->GetRadius()[0]*std::sin(angle));
                OutputImageType::RegionType outputRegion =
                        localOutputImage->GetLargestPossibleRegion();
                if( outputRegion.IsInside( localIndex ) )
                {
                    localOutputImage->SetPixel( localIndex, 255 );
                }
            }
            itCircles++;
        }

        /*try
         {
         writer->SetFileName(volume_out);
         writer->SetInput(localOutputImage);
         writer->Update();
         }
       catch( itk::ExceptionObject & excep )
         {
         std::cerr << "Exception caught !" << std::endl;
         std::cerr << excep << std::endl;
         }
*/
        typedef  itk::ImageFileWriter< ImageType  > WriterType;
        WriterType::Pointer writerHough = WriterType::New();

        typedef itk::Image<unsigned char, 2>  UnsignedCharImageType;
        //typedef itk::Image<float, 2>  FloatImageType;
        typedef  itk::ImageFileWriter<UnsignedCharImageType> WriterTypeCine;
        WriterTypeCine::Pointer writerCine = WriterTypeCine::New();

        stringstream stringFileHough;

        string typeTiff = ".tif";

        stringFileHough<<pathHough<<(i+1)<<typeTiff;

        string houghFile = stringFileHough.str();
        stringFileHough.str("");

        stringstream stringFileCine;

        stringFileHough<<pathCine<<(i+1)<<typeTiff;

        string cineFile = stringFileHough.str();
        stringFileHough.str("");

        typedef itk::RescaleIntensityImageFilter< ImageType2D, ImageType2D > RescaleFilterType2;
        RescaleFilterType2::Pointer rescaleFilter2 = RescaleFilterType2::New();
        rescaleFilter2->SetInput(extractImg->GetOutput());
        rescaleFilter2->SetOutputMinimum(0);
        rescaleFilter2->SetOutputMaximum(255);
        rescaleFilter2->Update();
        try{
            writerHough->SetFileName(houghFile);
            writerHough->SetInput(localOutputImage );
            writerHough->Update();

            writerCine->SetFileName(cineFile);
            writerCine->SetInput(rescaleFilter2->GetOutput());
            writerCine->Update();
        }
        catch( itk::ExceptionObject & excep )
        {
            std::cerr << "Exception caught !" << std::endl;
            std::cerr << excep << std::endl;
        }
    }
    endocardiumFile.close();
    radiusFile.close();

    /*typedef itk::Image<unsigned char,3> ImageType2;
      ImageType2::Pointer imag_out = ImageType::New();
      imag_out = localOutputImage;
      return imag_out;*/

}
