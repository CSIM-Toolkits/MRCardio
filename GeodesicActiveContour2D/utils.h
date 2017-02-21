#ifndef UTILS_H
#define UTILS_H

#include "itkImage.h"

class Utils
{
public:
    Utils();

    typedef float PixelTypeFloat;
    typedef unsigned char PixelTypeUC;
    typedef itk::Image<PixelTypeFloat,2> ImageType;

    double GetPixel(ImageType::Pointer image, double x, double y);
    double GetPixel(ImageType::Pointer image, double x, double y, double z);
    double GetMinimum(ImageType::Pointer image);
    double GetMaximum(ImageType::Pointer image);
    double GetStd(ImageType::Pointer image);
    double GetMean(ImageType::Pointer image);
    double GetHeight(ImageType::Pointer image);
    double GetWidth(ImageType::Pointer image);
    double GetDepth(ImageType::Pointer image);
    void GetSeed(ImageType::Pointer image,int centerX, int centerY, int *x, int *y);
    void GetCenter(ImageType::Pointer image, int *x, int *y);
};

#endif // UTILS_H