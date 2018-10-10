#ifndef KMEANS_H
#define KMEANS_H

#include "itkImage.h"

using namespace std;

class kmeans
{
public:
    kmeans();
    typedef unsigned char PixelType;
    typedef itk::Image<PixelType,2> ImageType;
    void Execute(int first, int last);
    int t_final;

    struct passwd *pw;

    const char *homedir;
    string extractValues;
    string slices;
    string segmentedFinal;
    string pathSlices;
    string pathSegmentedFinal;
    string pathExtractValues;
    string pathKmeans;
    string km;


};

#endif // KMEANS_H
