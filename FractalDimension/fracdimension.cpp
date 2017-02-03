#include "fracdimension.h"

#include "itkImage.h"
#include <fstream>
//#include "itkImageToHistogramGenerator.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionConstIterator.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkFlatStructuringElement.h"
#include "itkSubtractImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

#include <iostream>
#include <string>
#include <math.h>
#include <fstream>

#include "itkBinaryBallStructuringElement.h"
#include <math.h>

using namespace std;

fracdimension::fracdimension()
{
}

/**
 * @brief fracdimension::GetBoxCountingDimension2D
 * Computing Box Counting for 2D images
 * @param image
 * @return dim
 */
double fracdimension::GetBoxCountingDimension2D(ImageType2D::Pointer image){

    typedef unsigned int PixelType;

    typedef itk::Image< PixelType, 2> ImageType;

    typedef itk::ShapedNeighborhoodIterator<ImageType> IteratorType;

    double dim = 0.0;
    int cont_A = 0;
    int cont_elem = 0;
    double sum_la = 0.0;
    bool isZero = false;
    bool isElement = false;
    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();
    double M = region[0];
    double r = 0;
    int k = 1;
    double *vetNR = new double[10];
    double *vetR = new double[10];

    double s = (pow((1+(k*2)),2)/2);
    while(s < (M)){
        cont_A = 0;
        dim = 0;
        int cont = 0;
        for(int a=0; a<region[0]; a = a + (1+((k-1)*2))){
            for(int b= 0; b<region[1]; b = b + (1+((k-1)*2))){
                isElement = false;
                for(int c = 0; c<(1+((k-1)*2)); c++){
                    for(int d = 0; d<(1+((k-1)*2)); d++){
                        const ImageType::IndexType index = {{a+c,b+d}};
                        if(((a+c) >= 0) && ((a+c) <= region[0]) && ((b+d) >= 0) && ((b+d) <= region[1])){
                            if(image->GetPixel(index) > 0){
                                isElement = true;
                            }
                        }
                    }
                }
                if(isElement){
                    cont_A++;
                }
                cont++;
            }
        }

        r = (1+((k-1)*2))/M;
        vetNR[k] = (log(cont_A));
        vetR[k] = (log(r));
        k++;
        s = (pow((1+(k*2)),2)/2);

    }
    double m,b;
    linreg(k-1,vetR,vetNR,&m,&b);
    free(vetR);
    free(vetNR);
    dim = -m;
    return dim;
}

/**
 * @brief fracdimension::GetBoxCountingDimension3D
 * Computing Box Counting for 3D images
 * @param image
 * @return dim
 */
double fracdimension::GetBoxCountingDimension3D(ImageType::Pointer image){
    typedef unsigned int PixelType;

    typedef itk::Image< PixelType, 3> ImageType;

    typedef itk::ShapedNeighborhoodIterator<ImageType> IteratorType;

    double dim = 0.0;
    int cont_A = 0;
    int cont_elem = 0;
    double sum_la = 0.0;
    bool isZero = false;
    bool isElement = false;
    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();
    double M = region[0];
    double r = 0;
    int k = 1;
    double *vetNR = new double[10];
    double *vetR = new double[10];

    double s = (pow((1+(k*2)),2)/2);
    while(s < (M)){
        cont_A = 0;
        dim = 0;
        int cont = 0;
        for(int a=0; a<region[0]; a = a + (1+((k-1)*2))){
            for(int b= 0; b<region[1]; b = b + (1+((k-1)*2))){
                for(int u= 0; u<region[2]; u = u + (1+((k-1)*2))){
                    isElement = false;
                    for(int c = 0; c<(1+((k-1)*2)); c++){
                        for(int d = 0; d<(1+((k-1)*2)); d++){
                            for(int e = 0; e<(1+((k-1)*2)); e++){
                                const ImageType::IndexType index = {{a+c,b+d,u+e}};
                                if(((a+c) >= 0) && ((a+c) <= region[0]) && ((b+d) >= 0) && ((b+d) <= region[1]) && ((u+e) >= 0) && ((u+e) <= region[2])){
                                    if(image->GetPixel(index) > 0){
                                        isElement = true;
                                    }
                                }
                            }
                        }
                    }
                    if(isElement){
                        cont_A++;
                    }
                    cont++;
                }
            }
        }

        r = (1+((k-1)*2))/M;
        vetNR[k] = (log(cont_A));
        vetR[k] = (log(r));
        k++;
        s = (pow((1+(k*2)),2)/2);

    }
    double m,b;
    linreg(k-1,vetR,vetNR,&m,&b);
    free(vetR);
    free(vetNR);
    dim = -m;
    return dim;
}

double fracdimension::GetDBCDimension(ImageType::Pointer Image){
    double dim = 0.0;
    int r = 2;
    int cont_A = 0;
    int cont_elem = 0;
    double sum_la = 0.0;
    int nr = 0;
    int min = 9999;
    int auxMin = 0;
    int auxMax = 0;
    int max = 0;
    typedef unsigned int PixelType;
    typedef itk::Image<PixelType,3> ImageType;
    typedef itk::ShapedNeighborhoodIterator<ImageType> IteratorType;
    typedef itk::ConstNeighborhoodIterator<ImageType> ShapedNeighborhoodIteratorType;
    typedef itk::FlatStructuringElement<3> FlatStructuringElementType;
    FlatStructuringElementType::RadiusType r1;
    r1.Fill(r);
    //ShapedNeighborhoodIteratorType::RadiusType radius;
    //radius.Fill(3);
    itk::Size<3> radius;
    radius.Fill(2);

    IteratorType iterator(radius,Image,Image->GetLargestPossibleRegion());

    IteratorType::OffsetType top = {{0,-1}};
    iterator.ActivateOffset((top));
    //IteratorType::OffsetType bottom = {{0,1}};
    //iterator.ActivateOffset(bottom);
    IteratorType::OffsetType left = {{-1,0}};
    iterator.ActivateOffset(left);
    //IteratorType::OffsetType right = {{1,0}};
    //iterator.ActivateOffset((right));

    //IteratorType::OffsetType topr = {{1,-1}};
    //iterator.ActivateOffset((topr));
    //IteratorType::OffsetType bottomr = {{1,1}};
    //iterator.ActivateOffset(bottomr);
    IteratorType::OffsetType topl = {{-1,-1}};
    iterator.ActivateOffset(topl);
    //IteratorType::OffsetType bottoml = {{-1,1}};
    //iterator.ActivateOffset((bottoml));
    IteratorType::OffsetType center = {{0,0}};
    iterator.ActivateOffset((center));

    /*IteratorType::OffsetType one = {{0,-2}};
    iterator.ActivateOffset((one));
    IteratorType::OffsetType two = {{2,0}};
    iterator.ActivateOffset((two));
    IteratorType::OffsetType three = {{0,2}};
    iterator.ActivateOffset((three));
    IteratorType::OffsetType four = {{-2,0}};
    iterator.ActivateOffset((four));*/
    ofstream myfile ("/home/gustavo/temp/dimensions.txt");
    if (myfile.is_open())
    {

        //IteratorType::ConstIterator es = iterator.Begin();
        for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator){
            int cont_elem = 0;
            IteratorType::ConstIterator es = iterator.Begin();
            while(!es.IsAtEnd()){

                auxMin = es.Get();
                auxMax = es.Get();
                if(auxMin < min){
                    min = auxMin;
                }
                if(auxMax > max){
                    max = auxMax;
                }
                cont_elem++;
                es++;
                //cout<<"MIN: "<<min<<"  Element: "<<cont_elem<<endl;
                //cout<<"MAX: "<<max<<endl;

                //cout<<"CONT: "<<cont_elem<<endl;
            }
            nr = nr + (((int)(max/2) - (int)(min/2)) + 1);
            myfile<<"NR: "<<nr<<endl;
            //cout<<"NR: "<<auxMin<<endl;
            min = 9999;
            auxMin = 0;
            auxMax = 0;
            max = 0;

        }
    }

    dim = ((log(nr))/(2.10));
    myfile.close();

    return dim;

}

void fracdimension::linreg(int n, double x[], double y[], double* m, double* b)
{
    double   sumx = 0.0;                        /* sum of x                      */
    double   sumx2 = 0.0;                       /* sum of x**2                   */
    double   sumxy = 0.0;                       /* sum of x * y                  */
    double   sumy = 0.0;                        /* sum of y                      */
    double   sumy2 = 0.0;                       /* sum of y**2                   */

    for (int i=1;i<=n;i++)
    {
        sumx  += x[i];
        sumx2 += pow(x[i],2);
        sumxy += x[i] * y[i];
        sumy  += y[i];
        sumy2 += pow(y[i],2);
    }

    double denom = (n * sumx2 - pow(sumx,2));
    if (denom == 0) {
        // singular matrix. can't solve the problem.
        *m = 0;
        *b = 0;
    }

    *m = (n * sumxy  -  sumx * sumy) / denom;
    *b = (sumy * sumx2  -  sumx * sumxy) / denom;

}
