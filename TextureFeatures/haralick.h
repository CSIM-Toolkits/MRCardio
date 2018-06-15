#ifndef HARALICK_H
#define HARALICK_H

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkHistogramToTextureFeaturesFilter.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkVectorContainer.h"
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"

class Haralick
{
public:
    Haralick();

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

    typedef float PixelType;
    typedef itk::Image<PixelType,3> ImageType;
    void calcTextureFeatureImage(OffsetType offset,
                                 InternalImageType::Pointer inputImage, InternalImageType::Pointer outInertia,
                                 InternalImageType::Pointer outCorrelation, InternalImageType::Pointer outEnergy);
    void Execute(ImageType::Pointer Image);
};

#endif // HARALICK_H
