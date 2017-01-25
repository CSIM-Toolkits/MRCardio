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
    glcmGenerator3D->SetNumberOfBinsPerAxis(16); //reasonable number of bins
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


