#ifndef FRACDIMENSION_H
#define FRACDIMENSION_H

#include <itkImage.h>
#include <math.h>
#define REAL double
inline static REAL sqr(REAL x) {
        return x*x;
    }

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
    int linreg(int n, const REAL x[], const REAL y[], REAL* m, REAL* b, REAL* r);
};

#endif // FRACDIMENSION_H
