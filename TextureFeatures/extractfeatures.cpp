#include "extractfeatures.h"
#include <itkImageFileReader.h>
#include <itkLabelImageToStatisticsLabelMapFilter.h>

#include <itkImageFileWriter.h>
#include "itkImage.h"
#include <itkPasteImageFilter.h>
#include <itkShapeLabelMapFilter.hxx>
#include <itkGeometryUtilities.h>
#include <itkScalarImageToTextureFeaturesFilter.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <itkDenseFrequencyContainer2.h>
#include "itkHistogramToTextureFeaturesFilter.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkVectorContainer.h"
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"

ExtractFeatures::ExtractFeatures()
{
}

typedef itk::Image<float, 3> InternalImageType;
typedef itk::Image<unsigned char, 3> VisualizingImageType;
typedef itk::Neighborhood<float, 3> NeighborhoodType;
typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<InternalImageType>
Image2CoOccuranceType;
typedef Image2CoOccuranceType::HistogramType HistogramType;
typedef itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType> Hist2FeaturesType;
typedef InternalImageType::OffsetType OffsetType;
typedef itk::AddImageFilter <InternalImageType> AddImageFilterType;
typedef itk::MultiplyImageFilter<InternalImageType> MultiplyImageFilterType;
int ExtractFeatures::Extract(OffsetType offset, InternalImageType::Pointer inputImage)
{    
    // principal variables
    //Gray Level Co-occurance Matrix Generator
    Image2CoOccuranceType::Pointer glcmGenerator=Image2CoOccuranceType::New();
    glcmGenerator->SetOffset(offset);
    glcmGenerator->SetNumberOfBinsPerAxis(16); //reasonable number of bins
    glcmGenerator->SetPixelValueMinMax(0, 255); //for input UCHAR pixel type
    Hist2FeaturesType::Pointer featureCalc=Hist2FeaturesType::New();
    //Region Of Interest
    typedef itk::RegionOfInterestImageFilter<InternalImageType,InternalImageType> roiType;
    roiType::Pointer roi=roiType::New();
    roi->SetInput(inputImage);

    InternalImageType::RegionType window;
    InternalImageType::RegionType::SizeType size;
    size.Fill(50);
    window.SetSize(size);

    window.SetIndex(0,0);
    window.SetIndex(1,0);
    window.SetIndex(2,0);

    roi->SetRegionOfInterest(window);
    roi->Update();

    glcmGenerator->SetInput(roi->GetOutput());
    glcmGenerator->Update();

    featureCalc->SetInput(glcmGenerator->GetOutput());
    featureCalc->Update();

    std::cout<<"\n Entropy : ";
    std::cout<<featureCalc->GetEntropy()<<"\n Energy";
    std::cout<<featureCalc->GetEnergy()<<"\n Correlation";
    std::cout<<featureCalc->GetCorrelation()<<"\n Inertia";
    std::cout<<featureCalc->GetInertia()<<"\n HaralickCorrelation";
    std::cout<<featureCalc->GetHaralickCorrelation()<<"\n InverseDifferenceMoment";
    std::cout<<featureCalc->GetInverseDifferenceMoment()<<"\nClusterProminence";
    std::cout<<featureCalc->GetClusterProminence()<<"\nClusterShade";
    std::cout<<featureCalc->GetClusterShade();
}




