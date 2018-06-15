#include "utils.h"

using namespace std;

/**
 * @brief Utils::Utils
 */
Utils::Utils()
{

}

/**
 * @brief Utils::GetPixel
 * @param image
 * @param x
 * @param y
 * @return double
 */
double Utils::GetPixel(ImageType::Pointer image, double x, double y){
    ImageType::IndexType pixelIndex;

    pixelIndex[0] = x;      // x position of the pixel
    pixelIndex[1] = y;      // y position of the pixel
    pixelIndex[2] = 1;
    ImageType::PixelType pixelValue = image->GetPixel( pixelIndex );

    return pixelValue;
}

/**
 * @brief Utils::GetPixel
 * @param image
 * @param x
 * @param y
 * @param z
 * @return double
 */
double Utils::GetPixel(ImageType3D::Pointer image, double x, double y, double z){
    ImageType3D::IndexType pixelIndex;

    pixelIndex[0] = x;      // x position of the pixel
    pixelIndex[1] = y;      // y position of the pixel
    pixelIndex[2] = z;
    ImageType3D::PixelType pixelValue = image->GetPixel( pixelIndex );

    return pixelValue;
}

/**
 * @brief Utils::GetStd
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
 * @param image
 * @return double
 */
double Utils::GetHeight(ImageType::Pointer image){
    typedef itk::Image<unsigned char, 2>  ImageType;

    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[0];
}

/**
 * @brief Utils::GetWidth
 * @param image
 * @return double
 */
double Utils::GetWidth(ImageType::Pointer image){
    typedef itk::Image<unsigned char, 2>  ImageType;

    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[1];
}

/**
 * @brief Utils::GetHeight3D
 * @param image
 * @return double
 */
double Utils::GetHeight3D(ImageType3D::Pointer image){
    typedef itk::Image<unsigned char, 3>  ImageType3D;

    const ImageType3D::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[0];
}

/**
 * @brief Utils::GetWidth3D
 * @param image
 * @return double
 */
double Utils::GetWidth3D(ImageType3D::Pointer image){
    typedef itk::Image<unsigned char, 3>  ImageType3D;

    const ImageType3D::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[1];
}

/**
 * @brief Utils::GetDepth3D
 * @param image
 * @return double
 */
double Utils::GetDepth3D(ImageType3D::Pointer image){
    typedef itk::Image<unsigned char, 3>  ImageType3D;

    const ImageType3D::SizeType region = image->GetLargestPossibleRegion().GetSize();

    return region[2];
}

