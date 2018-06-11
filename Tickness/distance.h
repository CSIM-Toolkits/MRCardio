#ifndef DISTANCE_H
#define DISTANCE_H

#include <itkImage.h>
#include <math.h>
#include "QString"
using namespace std;

class Distance
{
public:
    Distance();
    typedef unsigned int PixelType;
    typedef itk::Image<PixelType,3> ImageType;
    typedef itk::Image<PixelType,2> ImageType2D;
    void GetTickness(int first, int last);
    double calcDistance(ImageType2D::Pointer image, int x, int y, int pos);

    struct passwd *pw;

    const char *homedir;

    string endocardium;
    string radius;
    string slices;
    string segmentedFinal;

    string pathEndocardium;
    string pathRadius;
    string pathSlices;
    string pathSegmentedFinal;

    string extractValues;
    string pathExtractValues;

};

#endif // DISTANCE_H
