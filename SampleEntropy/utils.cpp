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
    pixelIndex[2] = 1;
    ImageType::PixelType pixelValue = image->GetPixel( pixelIndex );

    return pixelValue;
}

double Utils::GetPixel(ImageType3D::Pointer image, double x, double y, double z){
    ImageType3D::IndexType pixelIndex;

    pixelIndex[0] = x;      // x position of the pixel
    pixelIndex[1] = y;      // y position of the pixel
    pixelIndex[2] = z;
    ImageType3D::PixelType pixelValue = image->GetPixel( pixelIndex );

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
double Utils::GetHeight3D(ImageType3D::Pointer image){
    typedef itk::Image<unsigned char, 3>  ImageType3D;

    const ImageType3D::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[0];
}

double Utils::GetWidth3D(ImageType3D::Pointer image){
    typedef itk::Image<unsigned char, 3>  ImageType3D;

    const ImageType3D::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[1];
}

double Utils::GetDepth3D(ImageType3D::Pointer image){
    typedef itk::Image<unsigned char, 3>  ImageType3D;

    const ImageType3D::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[2];
}

