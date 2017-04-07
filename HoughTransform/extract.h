#ifndef EXTRACT_H
#define EXTRACT_H

#include "itkImage.h"
#include "QString"
using namespace std;

class Extract
{
public:
    Extract();

    void Execute(int first,int last,const char* volume, int minimum, int maximum);

    struct passwd *pw;

    const char *homedir;

    string endocardium;
    string radius;
    string slices;
    string hough;
    string cine;

    string pathEndocardium;
    string pathRadius;
    string pathSlices;
    string pathHough;
    string pathCine;
};

#endif // EXTRACT_H
