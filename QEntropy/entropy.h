#ifndef ENTROPY_H
#define ENTROPY_H

#include "itkImage.h"

class entropy
{
public:
    entropy();
    typedef unsigned char PixelType;
    typedef itk::Image<PixelType,3> ImageType;
    double Execute(ImageType::Pointer Image, double q_entropy);
    int t_final;
};

#endif // ENTROPY_H
