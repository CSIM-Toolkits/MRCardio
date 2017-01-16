#ifndef EXTRACTFEATURES_H
#define EXTRACTFEATURES_H

#include <itkImageFileReader.h>
#include <itkLabelImageToStatisticsLabelMapFilter.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <itkDenseFrequencyContainer2.h>
#include "itkHistogramToTextureFeaturesFilter.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkVectorContainer.h"
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"

class ExtractFeatures
{
public:
    ExtractFeatures();
    //definitions of used types
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
    int Extract(OffsetType offset, InternalImageType::Pointer inputImage);
};

#endif // EXTRACTFEATURES_H
