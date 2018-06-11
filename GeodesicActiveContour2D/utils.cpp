#include "utils.h"
#include "itkImage.h"
#include "math.h"
#include "iostream"

using namespace std;

/**
 * @brief Utils::Utils
 */
Utils::Utils()
{

}

/**
 * @brief Utils::GetPixel
 * Return the pixel value in (x,y) coordinates
 * @param image
 * @param x
 * @param y
 * @return pixel
 */
double Utils::GetPixel(ImageType::Pointer image, double x, double y){
    ImageType::IndexType pixelIndex;

    pixelIndex[0] = x;      // x position of the pixel
    pixelIndex[1] = y;      // y position of the pixel
    pixelIndex[2] = 1;      // y position of the pixel
    ImageType::PixelType pixelValue = image->GetPixel( pixelIndex );

    return pixelValue;
}

/**
 * @brief Utils::GetPixel
 * Return the pixel value in (x,y,z) coordinates
 * @param image
 * @param x
 * @param y
 * @param z
 * @return pixel
 */
double Utils::GetPixel(ImageType::Pointer image, double x, double y, double z){
    ImageType::IndexType pixelIndex;

    pixelIndex[0] = x;      // x position of the pixel
    pixelIndex[1] = y;      // y position of the pixel
    pixelIndex[2] = z;
    ImageType::PixelType pixelValue = image->GetPixel( pixelIndex );

    return pixelValue;
}

/**
 * @brief Utils::GetStd
 * Return the standard deviation of image
 * @param image
 * @return double
 */
double Utils::GetStd(ImageType::Pointer image){
    double std = 0.0;
    double meanImage = GetMean(image);
    double width = GetWidth(image);
    double height = GetHeight(image);
    double M = (width * height) - 1;
    for(int i = 0; i < width; i++){
        for(int j = 0; j < height; j++){
            std += pow((GetPixel(image, i, j) - meanImage), 2);
        }
    }
    std /= M;
    return sqrt(std);
}

/**
 * @brief Utils::GetMean
 * Return the mean of image
 * @param image
 * @return double
 */
double Utils::GetMean(ImageType::Pointer image){
    double mean = 0.0;
    double width = GetWidth(image);
    double height = GetHeight(image);
    double N = (width * height) - 1;
    for(int i = 0; i < width; i++){
        for(int j = 0; j < height; j++){
            mean += GetPixel(image, i, j);
        }
    }

    return (mean/N);
}

/**
 * @brief Utils::GetHeight
 * Return the height of image
 * @param image
 * @return mixed
 */
double Utils::GetHeight(ImageType::Pointer image){
    typedef itk::Image<unsigned char, 2>  ImageType;

    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[0];
}

/**
 * @brief Utils::GetWidth
 * Return the width of image
 * @param image
 * @return mixed
 */
double Utils::GetWidth(ImageType::Pointer image){
    typedef itk::Image<unsigned char, 2>  ImageType;

    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[1];
}

/**
 * @brief Utils::GetDepth
 * Return the depth of image
 * @param image
 * @return mixed
 */
double Utils::GetDepth(ImageType::Pointer image){
    typedef itk::Image<unsigned char, 2>  ImageType;

    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[2];
}

/**
 * @brief Utils::GetSeedLeft
 * Find the seeds coordinates of septal myocardial position
 * @param image
 * @param centerX
 * @param centerY
 * @param x
 * @param y
 */
void Utils::GetSeedLeft(ImageType::Pointer image, int centerX, int centerY, int *x, int *y){
    double value = 9999;
    bool a = false;
    for(int i = centerX; i > 0; i--){
        value = (int) GetPixel(image, i, centerY);
        int aux = i;
        if(value > 0 && !a){
            *x = aux;
            *y = centerY;
            a = true;
            /*if(((GetPixel(image, centerX, centerY) - GetPixel(image, i, centerY)) / GetPixel(image, centerX, centerY)) > 0.4 ){
                a = true;
            }*/
        }
    }
}

/**
 * @brief Utils::GetSeedHight
 * Find the seeds coordinates of lateral myocardial position
 * @param image
 * @param centerX
 * @param centerY
 * @param x
 * @param y
 */
void Utils::GetSeedHight(ImageType::Pointer image, int centerX, int centerY, int *x, int *y){
    double height = GetHeight(image);
    double value = 9999;
    bool a = false;
    for(int i = centerX; i < height; i++){
        value = (int) GetPixel(image, i, centerY);
        int aux = i;
        if(value > 0 && !a){
            *x = aux;
            *y = centerY;
            a = true;
            /*if(((GetPixel(image, centerX, centerY) - GetPixel(image, i, centerY)) / GetPixel(image, centerX, centerY)) > 0.4 ){
                a = true;
            }*/
        }
    }
}

/**
 * @brief Utils::GetSeedUp
 * Find the seeds coordinates of anterior myocardial position
 * @param image
 * @param centerX
 * @param centerY
 * @param x
 * @param y
 */
void Utils::GetSeedUp(ImageType::Pointer image, int centerX, int centerY, int *x, int *y){
    double value = 9999;
    bool a = false;
    for(int i = centerY; i > 0; i--){
        value = (int) GetPixel(image, centerX, i);
        int aux = i;
        if(value > 0 && !a){
            *x = centerX;
            *y = aux;
            a = true;
            /*if(((GetPixel(image, centerX, centerY) - GetPixel(image, centerX, i)) / GetPixel(image, centerX, centerY)) > 0.4 ){
                a = true;
            }*/
        }
    }
}

/**
 * @brief Utils::GetSeedDown
 * Find the seeds coordinates of inferior myocardial position
 * @param image
 * @param centerX
 * @param centerY
 * @param x
 * @param y
 */
void Utils::GetSeedDown(ImageType::Pointer image, int centerX, int centerY, int *x, int *y){
    double width = GetWidth(image);
    double value = 9999;
    bool a = false;
    for(int i = centerY; i < width ; i++){
        value = (int) GetPixel(image, centerX, i);
        int aux = i;
        if(value > 0 && !a){
            *x = centerX;
            *y = aux;
            a = true;
            /*if(((GetPixel(image, centerX, centerY) - GetPixel(image, centerX, i)) / GetPixel(image, centerX, centerY)) > 0.4 ){
                a = true;
            }*/
        }
    }
}

/**
 * @brief Utils::GetCenter
 * Get the center of image
 * @param image
 * @param x
 * @param y
 */
void Utils::GetCenter(ImageType::Pointer image, int *x, int *y){
    int width = (int) GetWidth(image)/2;
    int height = (int) GetHeight(image)/2;
    *x = width;
    *y = height;
}

/**
 * @brief Utils::GetPerimeter
 * Get the perimeter of image contour
 * @param image
 * @return double
 */
double Utils::GetPerimeter(ImageTypeUC::Pointer image){
    double perimeter = 0.0;
    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();
    for(int c = 0; c<region[0]; c++){
        for(int d = 0; d<region[1]; d++){
            const ImageType::IndexType index = {{c,d}};
            if(image->GetPixel(index) > 0){
                perimeter++;
            }

        }
    }
    return perimeter;
}

/**
 * @brief Utils::GetArea
 * Get the area of image
 * @param image
 * @return double
 */
double Utils::GetArea(ImageTypeUC::Pointer image){
    double area;
    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();
    for(int c = 0; c<region[0]; c++){
        for(int d = 0; d<region[1]; d++){
            const ImageType::IndexType index = {{c,d}};
            if(image->GetPixel(index) > 0){
                area++;
            }

        }
    }
    return area;
}
