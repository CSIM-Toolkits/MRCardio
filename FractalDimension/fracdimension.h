#ifndef FRACDIMENSION_H
#define FRACDIMENSION_H

#include <itkImage.h>
#include <math.h>

class fracdimension
{
public:
    fracdimension();
    typedef unsigned int PixelType;
    typedef itk::Image<PixelType,3> ImageType;
    typedef itk::Image<PixelType,2> ImageType2D;
    double GetBoxCountingDimension2D(ImageType2D::Pointer image);
    double GetBoxCountingDimension3D(ImageType::Pointer image);
    double GetDBCDimension(ImageType::Pointer Image);
    void linreg(int n, double x[], double y[], double* m, double* b);
};

#endif // FRACDIMENSION_H
