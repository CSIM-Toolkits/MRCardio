#ifndef SEGMENT_H
#define SEGMENT_H
#include "itkImage.h"
#include "QString"
using namespace std;

class segment
{
public:
    segment();
    typedef float PixelType;
    typedef itk::Image<PixelType,3> ImageType;
    void InternalEC(int first,int last,double sigma,double sig_min, double sig_max, double propagation, double curvature, double advection, double rms, int iterations, double timestep, int it_dif, double conductance, double alpha, double meta, double distance);
    void MyocardiumEC(int first,int last,double sigma,double sig_min, double sig_max, double propagation, double curvature, double advection, double rms, int iterations, double timestep, int it_dif, double conductance, double alpha, double meta, double distance);

    void InternalELV(int first,int last,double sigma,double sig_min, double sig_max, double propagation, double curvature, double advection, double rms, int iterations, double timestep, int it_dif, double conductance, double alpha, double meta, double distance);
    void MyocardiumELV(int first,int last,double sigma,double sig_min, double sig_max, double propagation, double curvature, double advection, double rms, int iterations, double timestep, int it_dif, double conductance, double alpha, double meta, double distance);

    struct passwd *pw;

    const char *homedir;

    string endocardium;
    string radius;
    string slices;
    string segmented;
    string segmentedArea;
    string segmentedFinal;
    string cine;
    string gacFilter;
    string gacGradient;
    string gacSigmoid;
    string gacMap;
    string gacFilterFinal;
    string gacGradientFinal;
    string gacSigmoidFinal;
    string gacMapFinal;

    string gacGradientMHA;
    string gacSigmoidMHA;
    string gacMapMHA;
    string gacGradientMHAFinal;
    string gacSigmoidMHAFinal;
    string gacMapMHAFinal;


    string pathEndocardium;
    string pathRadius;
    string pathSlices;
    string pathSegmented;
    string pathSegmentedArea;
    string pathSegmentedFinal;
    string pathCine;
    string pathGacFilter;
    string pathGacGradient;
    string pathGacSigmoid;
    string pathGacMap;
    string pathGacFilterFinal;
    string pathGacGradientFinal;
    string pathGacSigmoidFinal;
    string pathGacMapFinal;

    string pathGacGradientMHA;
    string pathGacSigmoidMHA;
    string pathGacMapMHA;
    string pathGacGradientMHAFinal;
    string pathGacSigmoidMHAFinal;
    string pathGacMapMHAFinal;

    string extractValues;
    string extractValuesMyocardium;
    string pathExtractValues;
    string pathExtractValuesMyocardium;
};

#endif // SEGMENT_H
