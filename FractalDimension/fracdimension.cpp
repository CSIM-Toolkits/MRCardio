#include "fracdimension.h"

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

    double dim = 0.0;
    double cont_A = 0.0;
    bool isElement = false;
    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();
    double x = region[0];
    double y = region[1];
    double M = x;
    if(M < y){
        M = y;
    }
    double r = M;
    int k = 1;
    double *vetNR = new double[40];
    double *vetR = new double[40];

    while(r > 2){
        cont_A = 0.0;
        int cont = 0;
        for(int a = 0; a < int(x); a = a + int(r)){
            for(int b = 0; b < int(y); b = b + int(r)){
                isElement = false;
                for(int c = 0; c < r; c++){
                    for(int d = 0; d < r; d++){
                        const ImageType::IndexType index = {{a+c,b+d}};
                        if(((a+c) >= 0) && ((a+c) <= x) && ((b+d) >= 0) && ((b+d) <= y)){
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

        vetNR[k] = (log(cont_A));
        vetR[k] = (log(r));
        k++;
        r = r/2.0;

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

    double dim = 0.0;
    int cont_A = 0;
    bool isElement = false;
    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();
    double M = region[0];
    double r = 0;
    int k = 1;
    double *vetNR = new double[40];
    double *vetR = new double[40];

    double s = (1+((k-1)*2));
    while(s < (M)){
        cont_A = 0;
        dim = 0;
        int cont = 0;
        for(int a= 0; a<region[0]; a = a + (1+((k-1)*2))){
            for(int b = 0; b<region[1]; b = b + (1+((k-1)*2))){
                for(int u = 0; u<region[2]; u = u + (1+((k-1)*2))){
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
        s = (1+((k-1)*2));

    }
    double m,b;
    linreg(k-1,vetR,vetNR,&m,&b);
    free(vetR);
    free(vetNR);
    dim = -m;
    return dim;
}

/**
 * @brief fracdimension::GetDBCDimension2D
 * Computing DBC Dimension for 2D images
 * @param image
 * @return dimension
 */
double fracdimension::GetDBCDimension2D(ImageType2D::Pointer image){
    typedef unsigned int PixelType;

    typedef itk::Image< PixelType, 2> ImageType;
    double dim = 0.0;
    int cont_A = 0;
    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();
    double x = region[0];
    double y = region[1];
    double M = x;
    double r = 0.0;
    int k = 1;
    double *vetNR = new (nothrow) double[90];
    double *vetR = new (nothrow) double[90];
    double s = 2.0;
    double h = 0.0;
    while(s <= (M/2)){
        double min = 9999.0;
        double max = 0.0;
        cont_A = 0;
        int cont = 0;
        for(int a= 0; a < x; a = a + int(s)){
            for(int b = 0; b < y; b = b + int(s)){
                for(int c = 0; c < int(s); c++){
                    for(int d = 0; d < int(s); d++){
                        const ImageType::IndexType index = {{a+c,b+d}};
                        if(((a+c) > 0) && ((a+c) < x) && ((b+d) > 0) && ((b+d) < y)){
                            double minB = image->GetPixel(index);
                            double maxB = image->GetPixel(index);
                            if(minB < min){
                                min = minB;
                            }
                            if(maxB > max){
                                max = maxB;
                            }
                        }
                    }
                }
                if(max > 0){
                    h = (s*256)/M;
                    int positionMax = int(ceil(max/h));
                    int positionMin = int(ceil(min/h));
                    cont_A = cont_A + (((positionMax)-(positionMin)) +1);
                    cont++;
                }
            }
        }

        r = s/M;
        vetNR[k] = (log(cont_A));
        vetR[k] = (log(r));
        k++;
        s = s + 1.0;

    }
    double m,b;
    linreg(k-1,vetR,vetNR,&m,&b);
    delete[] vetR;
    delete [] vetNR;
    dim = -m;
    return dim;
}

/**
 * @brief fracdimension::GetMinkowskiDimension2D
 * Computing Minkowski Dimension for 2D images
 * @param image
 * @return dimension
 */
double fracdimension::GetMinkowskiDimension2D(ImageType2D::Pointer image){
    typedef unsigned int PixelType;

    typedef itk::Image< PixelType, 2> ImageType;

    double dim = 0.0;
    double cont_A = 0.0;
    double contElement = 0.0;
    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();
    double x = region[0];
    double y = region[1];
    double M = x;
    double r = 1.0;
    int k = 1;
    double *vetNR = new (nothrow) double[40];
    double *vetR = new (nothrow) double[40];

    while(r < (M/2.0)){
        cont_A = 0.0;
        contElement = 0.0;
        for(int a = 0; a < x; a++){
            for(int b = 0; b < y; b++){
                const ImageType::IndexType origin = {{a,b}};
                if(image->GetPixel(origin) > 0){
                    for(int c = int(-r); c <= int(r); c++){
                        for(int d = int(-r); d <= int(r); d++){
                            if(((a+c) >= 0) && ((a+c) <= x) && ((b+d) >= 0) && ((b+d) <= y)){
                                contElement++;
                            }
                        }
                    }
                }
                cont_A = cont_A + contElement;
            }
        }
        vetNR[k] = (log(cont_A));
        vetR[k] = (log(r));
        r = r * 2.0;
        k++;
    }

    double m,b;
    linreg(k-1,vetR,vetNR,&m,&b);
    delete[] vetR;
    delete [] vetNR;
    dim = 2-m;
    return dim;
}

/**
 * @brief fracdimension::linreg
 * Computing linear regression based on x[] and y[]
 * @param n
 * @param x
 * @param y
 * @param m
 * @param b
 */
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
    if (denom == 0.0) {
        // singular matrix. can't solve the problem.
        *m = 0.0;
        *b = 0.0;
    }

    *m = (n * sumxy  -  sumx * sumy) / denom;
    *b = (sumy * sumx2  -  sumx * sumxy) / denom;

}
