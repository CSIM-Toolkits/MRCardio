#ifndef FFT_H
#define FFT_H

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkWrapPadImageFilter.h"
#include "itkForwardFFTImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

class FFT
{
public:
    FFT();
    typedef float PixelType;
    typedef itk::Image<PixelType,3> ImageType;
    void Execute(ImageType::Pointer image);
};

#endif // FFT_H
