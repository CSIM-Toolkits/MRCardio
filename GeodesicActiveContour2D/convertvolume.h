#ifndef CONVERTVOLUME_H
#define CONVERTVOLUME_H

#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkTIFFImageIO.h"f"
#include "itkPNGImageIO.h"
#include <string>

class convertVolume
{
public:
    convertVolume();
    typedef unsigned char PixelType;
    typedef itk::Image<PixelType,3> ImageType;
    ImageType::Pointer Convert(int first, int last);
};

#endif // CONVERTVOLUME_H
