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
    int Extract(OffsetType offset, InternalImageType::Pointer inputImage, double *entropy, double *energy,
                double *correlation, double *inertia, double *haralickCorrelation, double *inverseDifferenceMoment,
                double *clusterProminence, double *clusterShade);

    typedef itk::Image<float, 3> InternalImageType3D;
    typedef itk::Image<unsigned char, 3> VisualizingImageType3D;
    typedef itk::Neighborhood<float, 3> NeighborhoodType3D;
    typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<InternalImageType3D>
    Image2CoOccuranceType3D;
    typedef Image2CoOccuranceType3D::HistogramType HistogramType3D;
    typedef itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType3D> Hist2FeaturesType3D;
    typedef InternalImageType3D::OffsetType OffsetType3D;
    typedef itk::AddImageFilter <InternalImageType3D> AddImageFilterType3D;
    typedef itk::MultiplyImageFilter<InternalImageType3D> MultiplyImageFilterType3D;
    int Extract3D(OffsetType3D offset3D, InternalImageType3D::Pointer inputImage3D, double *entropy, double *energy,
                  double *correlation, double *inertia, double *haralickCorrelation, double *inverseDifferenceMoment,
                  double *clusterProminence, double *clusterShade);
};

#endif // EXTRACTFEATURES_H
