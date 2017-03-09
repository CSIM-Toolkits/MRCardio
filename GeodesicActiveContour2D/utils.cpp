#include "utils.h"
#include "itkImage.h"
#include "math.h"
#include "iostream"

using namespace std;

Utils::Utils()
{

}

double Utils::GetPixel(ImageType::Pointer image, double x, double y){
    ImageType::IndexType pixelIndex;

    pixelIndex[0] = x;      // x position of the pixel
    pixelIndex[1] = y;      // y position of the pixel
    pixelIndex[2] = 1;      // y position of the pixel
    ImageType::PixelType pixelValue = image->GetPixel( pixelIndex );

    return pixelValue;
}

double Utils::GetPixel(ImageType::Pointer image, double x, double y, double z){
    ImageType::IndexType pixelIndex;

    pixelIndex[0] = x;      // x position of the pixel
    pixelIndex[1] = y;      // y position of the pixel
    pixelIndex[2] = z;
    ImageType::PixelType pixelValue = image->GetPixel( pixelIndex );

    return pixelValue;
}

double Utils::GetMaximum(ImageType::Pointer image){
    double maximum;

    return maximum;
}

double Utils::GetMinimum(ImageType::Pointer image){
    double minimum;

    return minimum;
}

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

double Utils::GetHeight(ImageType::Pointer image){
    typedef itk::Image<unsigned char, 2>  ImageType;

    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[0];
}

double Utils::GetWidth(ImageType::Pointer image){
    typedef itk::Image<unsigned char, 2>  ImageType;

    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[1];
}

double Utils::GetDepth(ImageType::Pointer image){
    typedef itk::Image<unsigned char, 2>  ImageType;

    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[2];
}

void Utils::GetSeedLeft(ImageType::Pointer image, int centerX, int centerY, int *x, int *y){
    double width = GetWidth(image);
    double height = GetHeight(image);
    double value = 9999;
    bool a = false;
    for(int i = centerX; i > 0; i--){
        value = (int) GetPixel(image, i, centerY);
        int aux = i;
        if(value > 0 && !a){
            *x = aux;
            *y = centerY;
            a = true;
        }
    }
}

void Utils::GetSeedHight(ImageType::Pointer image, int centerX, int centerY, int *x, int *y){
    double width = GetWidth(image);
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
        }
    }
}

void Utils::GetSeedUp(ImageType::Pointer image, int centerX, int centerY, int *x, int *y){
    double width = GetWidth(image);
    double height = GetHeight(image);
    double value = 9999;
    bool a = false;
    for(int i = centerY; i > 0; i--){
        value = (int) GetPixel(image, centerX, i);
        int aux = i;
        if(value > 0 && !a){
            *x = centerX;
            *y = aux;
            a = true;
        }
    }
}

void Utils::GetSeedDown(ImageType::Pointer image, int centerX, int centerY, int *x, int *y){
    double width = GetWidth(image);
    double height = GetHeight(image);
    double value = 9999;
    bool a = false;
    for(int i = centerY; i < width ; i++){
        value = (int) GetPixel(image, centerX, i);
        int aux = i;
        if(value > 0 && !a){
            *x = centerX;
            *y = aux;
            a = true;
        }
    }
}

void Utils::GetCenter(ImageType::Pointer image, int *x, int *y){
    int width = (int) GetWidth(image)/2;
    int height = (int) GetHeight(image)/2;
    *x = width;
    *y = height;
}

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
