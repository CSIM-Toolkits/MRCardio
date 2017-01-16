#ifndef EXTRACT_H
#define EXTRACT_H

#include "itkImage.h"
#include "QString"

class Extract
{
public:
    Extract();
    //typedef float PixelType;
    //typedef unsigned char PixelType2;
    //typedef itk::Image<PixelType,2> ImageType;
    //typedef itk::Image<PixelType2,2> ImageType2;
    //typedef itk::Image<unsigned char,3> ImageType;
    void Execute(int first,int last,const char* volume, const char* volume_out);
};

#endif // EXTRACT_H
