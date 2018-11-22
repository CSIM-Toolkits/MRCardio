#ifndef MAPPING_H
#define MAPPING_H

#include "itkImage.h"

using namespace std;

class Mapping
{
public:
    Mapping();
    typedef float PixelTypeFloat;
    typedef unsigned short PixelTypeUC;
    typedef itk::Image<PixelTypeUC,2> ImageType;
    typedef itk::Image<PixelTypeUC,3> ImageType3D;
    void calcMapping(ImageType::Pointer fixedImage, ImageType::Pointer movingImage, ImageType::Pointer segmentedImage, int index);

    struct passwd *pw;

    const char *homedir;

    string deformable;
    string vector;
    string field;
    string segmentedFinal;

    string pathDeformable;
    string pathVector;
    string pathField;
    string pathSegmentedFinal;
};

#endif // MAPPING_H
