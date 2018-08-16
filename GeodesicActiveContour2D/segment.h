#ifndef SEGMENT_H
#define SEGMENT_H

#include "itkImage.h"

using namespace std;

class segment
{
public:
    segment();
    typedef float PixelType;
    typedef itk::Image<PixelType,3> ImageType;

    void InternalEC(int first,int last,double sigma, float sig_min, float sig_max, float propagation,
        float curvature, float advection, double rms, unsigned long iterations, double timestep, unsigned long it_dif, double conductance,
        double alpha, double meta, const float distance);

    void MyocardiumEC(int first,int last,double sigma, float sig_min, float sig_max, float propagation,
        float curvature, float advection, double rms, unsigned long iterations, double timestep, unsigned long it_dif, double conductance,
        double alpha, double meta, const float distance, bool up, bool down, bool left, bool hight);

    void InternalELV(int first,int last,double sigma, float sig_min, float sig_max,
        float propagation, float curvature, float advection, double rms, unsigned long iterations,
        double timestep, unsigned long it_dif, double conductance, double alpha, double meta, const float distance);

    void MyocardiumELV(int first,int last,double sigma, float sig_min, float sig_max,
        float propagation, float curvature, float advection, double rms, unsigned long iterations,
        double timestep, unsigned long it_dif, double conductance, double alpha, double meta, const float distance,
        bool up, bool down, bool left, bool hight);

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
    string seeds;

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
    string pathSeeds;
};

#endif // SEGMENT_H
