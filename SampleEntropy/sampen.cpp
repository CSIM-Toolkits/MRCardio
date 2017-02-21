#include "sampen.h"

#include <iostream>
#include "QString"
#include "itkImage.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <stdlib.h>
#include "utils.h"

using namespace std;

SampEn::SampEn()
{

}

bool SampEn::isSimilar(ImageType::Pointer image, int x1, int y1, int x2, int y2, int m, double r){
    Utils utils;
    for(int y = 0; y < m; y++){
        for(int x = 0; x< m; x++){
            double diff = abs((utils.GetPixel(image, (x1 + x), (y1 + y))) - (utils.GetPixel(image, (x2 + x), (y2 + y))));
            if(diff >= r){
                return false;
            }
        }
    }
    return true;
}

bool SampEn::isSimilarNext(ImageType::Pointer image, int x1, int y1, int x2, int y2, int m, double r){
    Utils utils;
    double diff;
    for(int y = 0; y <= m; y++){
        diff = abs((utils.GetPixel(image, (x1 + m), (y1 + y))) - (utils.GetPixel(image, (x2 + m), (y2 + y))));
        if(diff >= r){
            return false;
        }
    }
    for(int x = 0; x <= m; x++){
        diff = abs((utils.GetPixel(image, (x1 + x), (y1 + m))) - (utils.GetPixel(image, (x2 + x), (y2 + m))));
        if(diff >= r){
            return false;
        }
    }
    return true;
}

bool SampEn::isSimilar3D(ImageType3D::Pointer image, int x1, int y1, int x2, int y2, int z1, int z2, int m, double r){
    Utils utils;
    for(int y = 0; y < m; y++){
        for(int x = 0; x < m; x++){
            for(int z = 0; z < m; z++){
                double diff = abs((utils.GetPixel(image, (x1 + x), (y1 + y), (z1 + z))) - (utils.GetPixel(image, (x2 + x), (y2 + y), (z2 + z))));
                if(diff >= r){
                    return false;
                }
            }
        }
    }
    return true;
}

bool SampEn::isSimilarNext3D(ImageType3D::Pointer image, int x1, int y1, int x2, int y2, int z1, int z2, int m, double r){
    Utils utils;
    double diff;
    for(int y = 0; y <= m; y++){
        diff = abs((utils.GetPixel(image, (x1 + m), (y1 + y), (z1 + m))) - (utils.GetPixel(image, (x2 + m), (y2 + y), (z2 + m))));
        if(diff >= r){
            return false;
        }
    }
    for(int x = 0; x <= m; x++){
        diff = abs((utils.GetPixel(image, (x1 + x), (y1 + m), (z1 + m))) - (utils.GetPixel(image, (x2 + x), (y2 + m), (z2 + m))));
        if(diff >= r){
            return false;
        }
    }
    for(int z = 0; z <= m; z++){
        diff = abs((utils.GetPixel(image, (x1 + m), (y1 + m), (z1 + z))) - (utils.GetPixel(image, (x2 + m), (y2 + m), (z2 + z))));
        if(diff >= r){
            return false;
        }
    }
    return true;
}

double SampEn::calcSampleEn2D(ImageType::Pointer image, int m, double r){
    Utils utils;
    double sampleEntropy;
    double tol = r;
    int nx = utils.GetWidth(image);
    int ny = utils.GetHeight(image);
    int A = 0;
    int B = 0;
    int Cim, Cim1;
    double Cm = 0.0;
    double Cm1 = 0.0;

    double den = (nx - m) * (ny - m);

    for(int yi = 0; yi < (ny - m); yi++){
        for(int xi = 0; xi < (nx - m); xi++){
            Cim = Cim1 = 0;
            int yj = yi;
            int xj = xi + 1;
            while (xj < (nx - m)) {
                if (isSimilar(image, xi, yi, xj, yj, m, tol)) {
                    B++;
                    Cim++;

                    if (isSimilarNext(image, xi, yi, xj, yj, m, tol)) {
                        A++;
                        Cim1++;
                    }
                }
                xj++;
            }

            for (yj = yi + 1; yj < (ny - m); yj++) {
                for (xj = 0; xj < (nx - m); xj++) {
                    if (isSimilar(image, xi, yi, xj, yj, m, tol)) {
                        B++;
                        Cim++;

                        if (isSimilarNext(image, xi, yi, xj, yj, m, tol)) {
                            A++;
                            Cim1++;
                        }
                    }
                }
            }

            Cm += Cim / (den - 1);
            Cm1 += Cim1 / (den - 1);
            //cout<<"CM: "<<Cm<<endl;
            //cout<<"CM1: "<<Cm1<<endl;
        }
    }
    Cm /= den;
    Cm1 /= den;

    sampleEntropy = -(log(((double) Cm1) / ((double) Cm)));
    cout<<"Results: "<<sampleEntropy<<endl;

    return sampleEntropy;
}

double SampEn::calcSampleEn3D(ImageType3D::Pointer image, int m, double r){
    Utils utils;
    double sampleEntropy;
    double tol = r;
    int nx = utils.GetWidth3D(image);
    int ny = utils.GetHeight3D(image);
    int nz = utils.GetDepth3D(image);
    int A = 0;
    int B = 0;
    int Cim, Cim1;
    double Cm = 0.0;
    double Cm1 = 0.0;    

    double den = (nx - m) * (ny - m) * (nz - m);

    for(int yi = 0; yi < (ny - m); yi++){
        for(int xi = 0; xi < (nx - m); xi++){
            for(int zi = 0; zi < (nz - m); zi++){
                Cim = Cim1 = 0;
                int zj = zi;
                int yj = yi;
                int xj = xi + 1;
                while (xj < (nx - m)) {
                    if (isSimilar3D(image, xi, yi, xj, yj, zi, zj, m, tol)) {
                        B++;
                        Cim++;

                        if (isSimilarNext3D(image, xi, yi, xj, yj, zi, zj, m, tol)) {
                            A++;
                            Cim1++;
                        }
                    }
                    xj++;
                }
                for (zj = zi + 1; zj < (nz - m); zj++) {
                    for (yj = yi + 1; yj < (ny - m); yj++) {
                        for (xj = 0; xj < (nx - m); xj++) {
                            if (isSimilar3D(image, xi, yi, xj, yj, zi, zj, m, tol)) {
                                B++;
                                Cim++;

                                if (isSimilarNext3D(image, xi, yi, xj, yj, zi, zj, m, tol)) {
                                    A++;
                                    Cim1++;
                                }
                            }
                        }
                    }
                }

                Cm += Cim / (den - 1);
                Cm1 += Cim1 / (den - 1);

            }
            cout<<"CM: "<<Cm<<endl;
            cout<<"CM1: "<<Cm1<<endl;
        }
    }
    Cm /= den;
    Cm1 /= den;

    sampleEntropy = -(log(((double) Cm1) / ((double) Cm)));
    cout<<"Results: "<<sampleEntropy<<endl;


    return sampleEntropy;
}
