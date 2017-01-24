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
    ofstream haralick2D("/home/gustavo/temp/haralick2D.txt");
    if (haralick2D.is_open())
    {
        haralick2D<<"Entropy: "<<featureCalc->GetEntropy()<<endl;
        haralick2D<<"Energy: "<<featureCalc->GetEnergy()<<endl;
        haralick2D<<"Correlation: "<<featureCalc->GetCorrelation()<<endl;
        haralick2D<<"Inertia: "<<featureCalc->GetInertia()<<endl;
        haralick2D<<"HaralickCorrelation: "<<featureCalc->GetHaralickCorrelation()<<endl;
        haralick2D<<"InverseDifferenceMoment: "<<featureCalc->GetInverseDifferenceMoment()<<endl;
        haralick2D<<"ClusterProminence: "<<featureCalc->GetClusterProminence()<<endl;
        haralick2D<<"ClusterShade: "<<featureCalc->GetClusterShade()<<endl;
    }
    haralick2D.close();
}

int ExtractFeatures::Extract3D(OffsetType3D offset3D, InternalImageType3D::Pointer inputImage3D)
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

    ofstream haralick3D("/home/gustavo/temp/haralick3D.txt");
    if (haralick3D.is_open())
    {
        haralick3D<<"Entropy: "<<featureCalc3D->GetEntropy()<<endl;
        haralick3D<<"Energy: "<<featureCalc3D->GetEnergy()<<endl;
        haralick3D<<"Correlation: "<<featureCalc3D->GetCorrelation()<<endl;
        haralick3D<<"Inertia: "<<featureCalc3D->GetInertia()<<endl;
        haralick3D<<"HaralickCorrelation: "<<featureCalc3D->GetHaralickCorrelation()<<endl;
        haralick3D<<"InverseDifferenceMoment: "<<featureCalc3D->GetInverseDifferenceMoment()<<endl;
        haralick3D<<"ClusterProminence: "<<featureCalc3D->GetClusterProminence()<<endl;
        haralick3D<<"ClusterShade: "<<featureCalc3D->GetClusterShade()<<endl;
    }
    haralick3D.close();
}


