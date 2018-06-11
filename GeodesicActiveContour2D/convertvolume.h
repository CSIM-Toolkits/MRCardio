#ifndef CONVERTVOLUME_H
#define CONVERTVOLUME_H

#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkTIFFImageIO.h"
#include "itkPNGImageIO.h"
#include <string>
using namespace std;

class convertVolume
{
public:
    convertVolume();
    typedef float PixelType;
    typedef itk::Image<PixelType,3> ImageType;
    ImageType::Pointer Convert(int first, int last);

    struct passwd *pw;

    const char *homedir;

    string segmentedFinal;
    string pathSegmentedFinal;
    string output;
    string pathOutput;

};

#endif // CONVERTVOLUME_H
