#ifndef FRACDIMENSION_H
#define FRACDIMENSION_H

#include <itkImage.h>
class fracdimension
{
public:
    fracdimension();
    typedef unsigned int PixelType;
    typedef itk::Image<PixelType,3> ImageType;
    double GetBoxCountingDimension2D(ImageType::Pointer image);
    double GetBoxCountingDimension3D(ImageType::Pointer image);
    void MorphologicalGradient(ImageType::Pointer Image);
    void GradientMagnitude(ImageType::Pointer Image, const char* volume_out);
    double GetDBCDimension(ImageType::Pointer Image);
    ImageType::Pointer edge;
};

#endif // FRACDIMENSION_H
