#ifndef FRACDIMENSION_H
#define FRACDIMENSION_H

#include <itkImage.h>
class fracdimension
{
public:
    fracdimension();
    typedef unsigned char PixelType;
    typedef itk::Image<PixelType,3> ImageType;
    double GetDimension(ImageType::Pointer Image);
    double GetMinkowskiDimension();
    void MorphologicalGradient(ImageType::Pointer Image);
    void GradientMagnitude(ImageType::Pointer Image, const char* volume_out);
    double GetDBCDimension(ImageType::Pointer Image);
    ImageType::Pointer edge;
};

#endif // FRACDIMENSION_H