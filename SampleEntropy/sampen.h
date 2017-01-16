#ifndef SAMPEN_H
#define SAMPEN_H

#include "itkImage.h"

class SampEn
{
public:
    SampEn();

    typedef float PixelTypeFloat;
    typedef unsigned char PixelTypeUC;
    typedef itk::Image<PixelTypeUC,3> ImageType;

    bool isSimilar(ImageType::Pointer image, int x1, int y1, int x2, int y2, int m, double r);
    bool isSimilarNext(ImageType::Pointer image, int x1, int y1, int x2, int y2, int m, double r);
    bool isSimilar3D(ImageType::Pointer image, int x1, int y1, int x2, int y2, int z1, int z2, int m, double r);
    bool isSimilarNext3D(ImageType::Pointer image, int x1, int y1, int x2, int y2, int z1, int z2, int m, double r);
    double calcSampleEn2D(ImageType::Pointer image, int m, double r);
    double calcSampleEn3D(ImageType::Pointer image, int m, double r);
};

#endif // SAMPEN_H
