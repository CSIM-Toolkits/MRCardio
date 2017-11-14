#ifndef CONVERT_H
#define CONVERT_H

#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkTIFFImageIO.h"
#include "itkPNGImageIO.h"
#include <string>
using namespace std;

class convert
{
public:
    convert();
    typedef float PixelType;
    typedef itk::Image<PixelType,3> ImageType;
    ImageType::Pointer ConvertImage(int first, int last);

    struct passwd *pw;

    const char *homedir;

    string segmentedFinal;
    string pathSegmentedFinal;
    string output;
    string pathOutput;

};

#endif // CONVERT_H
