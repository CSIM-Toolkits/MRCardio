#ifndef MAPPING_H
#define MAPPING_H

#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"

using namespace std;

class Mapping
{
public:
    Mapping();
    typedef float PixelTypeFloat;
    typedef unsigned short PixelTypeUC;
    typedef itk::Image<PixelTypeUC,2> ImageType;
    typedef itk::VectorImage<PixelTypeFloat,2> VectorImageType;
    typedef itk::Image<PixelTypeUC,3> ImageType3D;
    void calcMapping(ImageType::Pointer fixedImage, ImageType::Pointer movingImage, ImageType::Pointer segmentedImage, int index);
    void calcMagnitude(VectorImageType::Pointer vectorImage, int index);

    struct passwd *pw;

    const char *homedir;

    string deformable;
    string field;
    string segmentedFinal;
    string vector;
    string magnitude;

    string pathDeformable;
    string pathVector;
    string pathField;
    string pathSegmentedFinal;
    string PathVector;
    string PathMagnitude;
};

#endif // MAPPING_H
