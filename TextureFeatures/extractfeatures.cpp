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
#include <fstream>
#include <iostream>

using namespace std;

ExtractFeatures::ExtractFeatures()
{
}

typedef itk::Image<float, 2> InternalImageType;
typedef itk::Image<unsigned char, 2> VisualizingImageType;
typedef itk::Neighborhood<float, 2> NeighborhoodType;
typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<InternalImageType>
Image2CoOccuranceType;
typedef Image2CoOccuranceType::HistogramType HistogramType;
typedef itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType> Hist2FeaturesType;
typedef InternalImageType::OffsetType OffsetType;
typedef itk::AddImageFilter <InternalImageType> AddImageFilterType;
typedef itk::MultiplyImageFilter<InternalImageType> MultiplyImageFilterType;

typedef itk::Image<float, 3> InternalImageType3D;
typedef itk::Image<unsigned char, 3> VisualizingImageType3D;
typedef itk::Neighborhood<float, 3> NeighborhoodType3D;
typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<InternalImageType3D>
Image2CoOccuranceType3D;
typedef Image2CoOccuranceType::HistogramType HistogramType3D;
typedef itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType3D> Hist2FeaturesType3D;
typedef InternalImageType::OffsetType OffsetType3D;
typedef itk::AddImageFilter <InternalImageType3D> AddImageFilterType3D;
typedef itk::MultiplyImageFilter<InternalImageType3D> MultiplyImageFilterType3D;

int ExtractFeatures::Extract(OffsetType offset, InternalImageType::Pointer inputImage, double *entropy, double *energy,
                             double *correlation, double *inertia, double *haralickCorrelation, double *inverseDifferenceMoment,
                             double *clusterProminence, double *clusterShade)
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
    size.Fill(80);
    window.SetSize(size);

    window.SetIndex(0,0);
    window.SetIndex(1,0);

    roi->SetRegionOfInterest(window);
    roi->Update();

    glcmGenerator->SetInput(roi->GetOutput());
    glcmGenerator->Update();

    featureCalc->SetInput(glcmGenerator->GetOutput());
    featureCalc->Update();

    *entropy = featureCalc->GetEntropy();
    *energy = featureCalc->GetEnergy();
    *correlation = featureCalc->GetCorrelation();
    *inertia = featureCalc->GetInertia();
    *haralickCorrelation = featureCalc->GetHaralickCorrelation();
    *inverseDifferenceMoment = featureCalc->GetInverseDifferenceMoment();
    *clusterProminence = featureCalc->GetClusterProminence();
    *clusterShade = featureCalc->GetClusterShade();

}

int ExtractFeatures::Extract3D(OffsetType3D offset3D, InternalImageType3D::Pointer inputImage3D, double *entropy, double *energy,
                               double *correlation, double *inertia, double *haralickCorrelation, double *inverseDifferenceMoment,
                               double *clusterProminence, double *clusterShade)
{
    // principal variables
    //Gray Level Co-occurance Matrix Generator
    Image2CoOccuranceType3D::Pointer glcmGenerator3D=Image2CoOccuranceType3D::New();
    glcmGenerator3D->SetOffset(offset3D);
    glcmGenerator3D->SetNumberOfBinsPerAxis(256); //reasonable number of bins
    glcmGenerator3D->SetPixelValueMinMax(0, 255); //for input UCHAR pixel type
    Hist2FeaturesType3D::Pointer featureCalc3D=Hist2FeaturesType3D::New();
    //Region Of Interest
    typedef itk::RegionOfInterestImageFilter<InternalImageType3D,InternalImageType3D> roiType3D;
    roiType3D::Pointer roi3D=roiType3D::New();
    roi3D->SetInput(inputImage3D);

    InternalImageType3D::RegionType window3D;
    InternalImageType3D::RegionType::SizeType size3D;
    size3D.Fill(50);
    window3D.SetSize(size3D);

    window3D.SetIndex(0,0);
    window3D.SetIndex(1,0);
    window3D.SetIndex(2,0);

    roi3D->SetRegionOfInterest(window3D);
    roi3D->Update();

    glcmGenerator3D->SetInput(roi3D->GetOutput());
    glcmGenerator3D->Update();

    featureCalc3D->SetInput(glcmGenerator3D->GetOutput());
    featureCalc3D->Update();

    *entropy = featureCalc3D->GetEntropy();
    *energy = featureCalc3D->GetEnergy();
    *correlation = featureCalc3D->GetCorrelation();
    *inertia = featureCalc3D->GetInertia();
    *haralickCorrelation = featureCalc3D->GetHaralickCorrelation();
    *inverseDifferenceMoment = featureCalc3D->GetInverseDifferenceMoment();
    *clusterProminence = featureCalc3D->GetClusterProminence();
    *clusterShade = featureCalc3D->GetClusterShade();
}

int ExtractFeatures::ExtractRunLength2D(InternalImageType::Pointer inputImage2D, double *shortRunEmphasis, double *longRunEmphasis,
                                        double *greyLevelNonuniformity, double *runLengthNonuniformity, double *lowGrayLevelRunEmphasis,
                                        double *highGreyLevelRunEmphasis, double *shortRunLowGreyLevelEmphasis,
                                        double *shortRunHighGreyLevelEmphasis, double *longRunLowGreyLevelEmphasis,
                                        double *longRunHighGreyLevelEmphasis)
{
    const unsigned int ImageDimension = 2;
    typedef float PixelType;
    typedef float RealType;

    typedef itk::Image<PixelType, ImageDimension> ImageType;
    typedef itk::Image<RealType, ImageDimension> RealImageType;


    typedef itk::Statistics::DenseFrequencyContainer2 HistogramFrequencyContainerType;

    typedef itk::Statistics::ScalarImageToRunLengthFeaturesFilter
            <RealImageType, HistogramFrequencyContainerType> RunLengthFilterType;
    RunLengthFilterType::Pointer runLengthFilter = RunLengthFilterType::New();
    runLengthFilter->SetInput(inputImage2D);

    ImageType::Pointer mask = NULL;
    PixelType label = itk::NumericTraits<PixelType>::One;


    unsigned int numberOfBins = 256;
    runLengthFilter->SetNumberOfBinsPerAxis( numberOfBins );


    itk::ImageRegionIteratorWithIndex<ImageType> ItI( inputImage2D,
                                                      inputImage2D->GetLargestPossibleRegion() );

    PixelType maxValue = itk::NumericTraits<PixelType>::NonpositiveMin();
    PixelType minValue = itk::NumericTraits<PixelType>::max();

    typedef itk::BoundingBox<unsigned long,
            ImageDimension, RealType> BoundingBoxType;
    BoundingBoxType::Pointer bbox = BoundingBoxType::New();
    BoundingBoxType::PointsContainerPointer points
            = BoundingBoxType::PointsContainer::New();
    itk::Point<RealType, ImageDimension> point;

    unsigned int idx = 0;

    for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
        if ( !mask || ( mask->GetPixel( ItI.GetIndex() ) == label ) )
        {
            if ( ItI.Get() < minValue )
            {
                minValue = ItI.Get();
            }
            else if ( ItI.Get() > maxValue )
            {
                maxValue = ItI.Get();
            }
            inputImage2D->TransformIndexToPhysicalPoint( ItI.GetIndex(), point );
            points->InsertElement( idx++, point );
        }
    }
    bbox->SetPoints( points );
    bbox->ComputeBoundingBox();
    BoundingBoxType::PointType pointMin = bbox->GetMinimum();
    BoundingBoxType::PointType pointMax = bbox->GetMaximum();

    runLengthFilter->SetPixelValueMinMax( minValue, maxValue );
    runLengthFilter->SetDistanceValueMinMax( 0, pointMin.EuclideanDistanceTo( pointMax ) );
    runLengthFilter->SetNumberOfBinsPerAxis( numberOfBins );
    runLengthFilter->FastCalculationsOff();

    try
    {
        runLengthFilter->Update();

        RunLengthFilterType::FeatureValueVectorPointer means =
                runLengthFilter->GetFeatureMeans();
        const RunLengthFilterType::FeatureNameVector* names =
                runLengthFilter->GetRequestedFeatures();

        RunLengthFilterType::FeatureValueVector::ConstIterator mIt =
                means->Begin();
        RunLengthFilterType::FeatureNameVector::ConstIterator nIt =
                names->Begin();

        *shortRunEmphasis = mIt.Value(); ++mIt;
        *longRunEmphasis = mIt.Value(); ++mIt;
        *greyLevelNonuniformity = mIt.Value(); ++mIt;
        *runLengthNonuniformity = mIt.Value(); ++mIt;
        *lowGrayLevelRunEmphasis = mIt.Value(); ++mIt;
        *highGreyLevelRunEmphasis = mIt.Value(); ++mIt;
        *shortRunLowGreyLevelEmphasis = mIt.Value(); ++mIt;
        *shortRunHighGreyLevelEmphasis = mIt.Value(); ++mIt;
        *longRunLowGreyLevelEmphasis = mIt.Value(); ++mIt;
        *longRunHighGreyLevelEmphasis = mIt.Value(); ++mIt;

        return EXIT_SUCCESS;
    }
    catch(...)
    {
        return EXIT_FAILURE;
    }
}

int ExtractFeatures::ExtractRunLength3D(InternalImageType3D::Pointer inputImage3D, double *shortRunEmphasis, double *longRunEmphasis,
                                        double *greyLevelNonuniformity, double *runLengthNonuniformity, double *lowGrayLevelRunEmphasis,
                                        double *highGreyLevelRunEmphasis, double *shortRunLowGreyLevelEmphasis,
                                        double *shortRunHighGreyLevelEmphasis, double *longRunLowGreyLevelEmphasis,
                                        double *longRunHighGreyLevelEmphasis)
{
    const unsigned int ImageDimension = 3;
    typedef float PixelType;
    typedef float RealType;

    typedef itk::Image<PixelType, ImageDimension> ImageType;
    typedef itk::Image<RealType, ImageDimension> RealImageType;


    typedef itk::Statistics::DenseFrequencyContainer2 HistogramFrequencyContainerType;

    typedef itk::Statistics::ScalarImageToRunLengthFeaturesFilter
            <RealImageType, HistogramFrequencyContainerType> RunLengthFilterType;
    RunLengthFilterType::Pointer runLengthFilter = RunLengthFilterType::New();
    runLengthFilter->SetInput(inputImage3D);

    ImageType::Pointer mask = NULL;
    PixelType label = itk::NumericTraits<PixelType>::One;


    unsigned int numberOfBins = 256;
    runLengthFilter->SetNumberOfBinsPerAxis( numberOfBins );


    itk::ImageRegionIteratorWithIndex<ImageType> ItI( inputImage3D,
                                                      inputImage3D->GetLargestPossibleRegion() );

    PixelType maxValue = itk::NumericTraits<PixelType>::NonpositiveMin();
    PixelType minValue = itk::NumericTraits<PixelType>::max();

    typedef itk::BoundingBox<unsigned long,
            ImageDimension, RealType> BoundingBoxType;
    BoundingBoxType::Pointer bbox = BoundingBoxType::New();
    BoundingBoxType::PointsContainerPointer points
            = BoundingBoxType::PointsContainer::New();
    itk::Point<RealType, ImageDimension> point;

    unsigned int idx = 0;

    for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
        if ( !mask || ( mask->GetPixel( ItI.GetIndex() ) == label ) )
        {
            if ( ItI.Get() < minValue )
            {
                minValue = ItI.Get();
            }
            else if ( ItI.Get() > maxValue )
            {
                maxValue = ItI.Get();
            }
            inputImage3D->TransformIndexToPhysicalPoint( ItI.GetIndex(), point );
            points->InsertElement( idx++, point );
        }
    }
    bbox->SetPoints( points );
    bbox->ComputeBoundingBox();
    BoundingBoxType::PointType pointMin = bbox->GetMinimum();
    BoundingBoxType::PointType pointMax = bbox->GetMaximum();

    runLengthFilter->SetPixelValueMinMax( minValue, maxValue );
    runLengthFilter->SetDistanceValueMinMax( 0, pointMin.EuclideanDistanceTo( pointMax ) );
    runLengthFilter->SetNumberOfBinsPerAxis( numberOfBins );
    runLengthFilter->FastCalculationsOff();

    try
    {
        runLengthFilter->Update();

        RunLengthFilterType::FeatureValueVectorPointer means =
                runLengthFilter->GetFeatureMeans();
        const RunLengthFilterType::FeatureNameVector* names =
                runLengthFilter->GetRequestedFeatures();

        RunLengthFilterType::FeatureValueVector::ConstIterator mIt =
                means->Begin();
        RunLengthFilterType::FeatureNameVector::ConstIterator nIt =
                names->Begin();

        *shortRunEmphasis = mIt.Value(); ++mIt;
        *longRunEmphasis = mIt.Value(); ++mIt;
        *greyLevelNonuniformity = mIt.Value(); ++mIt;
        *runLengthNonuniformity = mIt.Value(); ++mIt;
        *lowGrayLevelRunEmphasis = mIt.Value(); ++mIt;
        *highGreyLevelRunEmphasis = mIt.Value(); ++mIt;
        *shortRunLowGreyLevelEmphasis = mIt.Value(); ++mIt;
        *shortRunHighGreyLevelEmphasis = mIt.Value(); ++mIt;
        *longRunLowGreyLevelEmphasis = mIt.Value(); ++mIt;
        *longRunHighGreyLevelEmphasis = mIt.Value(); ++mIt;

        return EXIT_SUCCESS;
    }
    catch(...)
    {
        return EXIT_FAILURE;
    }
}

