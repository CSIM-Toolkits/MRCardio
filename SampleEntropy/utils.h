#ifndef UTILS_H
#define UTILS_H

#include "itkImage.h"

class Utils
{
public:
    Utils();

    typedef float PixelTypeFloat;
    typedef unsigned char PixelTypeUC;
    typedef itk::Image<PixelTypeUC,2> ImageType;
    typedef itk::Image<PixelTypeUC,3> ImageType3D;

    double GetPixel(ImageType::Pointer image, double x, double y);
    double GetPixel(ImageType3D::Pointer image, double x, double y, double z);
    double GetStd(ImageType::Pointer image);
    double GetMean(ImageType::Pointer image);
    double GetHeight(ImageType::Pointer image);
    double GetWidth(ImageType::Pointer image);
    double GetHeight3D(ImageType3D::Pointer image);
    double GetWidth3D(ImageType3D::Pointer image);
    double GetDepth3D(ImageType3D::Pointer image);
};

#endif // UTILS_H
