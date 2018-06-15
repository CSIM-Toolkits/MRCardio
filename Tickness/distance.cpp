#include "distance.h"
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
#include "itkBinaryBallStructuringElement.h"
#include <pwd.h>

using namespace std;

/**
 * @brief Distance::Distance
 */
Distance::Distance()
{
    this->pw = getpwuid(getuid());

    this->homedir = this->pw->pw_dir;

    this->endocardium = "/temp/endocardium.txt";
    this->radius = "/temp/radius.txt";
    this->slices = "/temp/slices.txt";
    this->segmentedFinal = "/temp/segmentedFinal_";
    this->extractValues = "/temp/extractedValuesInternal.txt";

    this->pathExtractValues = this->homedir + this->extractValues;
    this->pathEndocardium = this->homedir + this->endocardium;
    this->pathSegmentedFinal = this->homedir + this->segmentedFinal;
}

/**
 * @brief Distance::GetTickness
 * @param first
 * @param last
 */
void Distance::GetTickness(int first, int last){

    string typeTiff = ".tif";
    ifstream endocardiumFile (this->pathEndocardium.c_str());
    ofstream extractValues(this->pathExtractValues.c_str());
    if (extractValues.is_open()){
        for(int i =first; i<last;i++){

            typedef   unsigned int           InternalPixelType;
            const     unsigned int    Dimension = 2;

            typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;

            typedef  itk::ImageFileReader< InternalImageType > ReaderType;
            ReaderType::Pointer reader = ReaderType::New();

            string line;
            int x = 0;
            int y = 0;

            if (endocardiumFile.is_open())
            {
                getline (endocardiumFile,line);
                string cmd = line;
                string arg;
                string::size_type pos = cmd.find(' ');
                if(cmd.npos != pos) {
                    arg = cmd.substr(pos + 1);
                    cmd = cmd.substr(0, pos);
                }
                x = atoi(cmd.c_str());
                y = atoi(arg.c_str());
            }

            stringstream stringFileSegmented;

            if(i<9)
                stringFileSegmented<<this->pathSegmentedFinal<<"00"<<(i+1)<<typeTiff;
            if(i>=9 && i<99)
                stringFileSegmented<<this->pathSegmentedFinal<<"0"<<(i+1)<<typeTiff;
            if(i>=99)
                stringFileSegmented<<this->pathSegmentedFinal<<(i+1)<<typeTiff;
            string segmentedFile = stringFileSegmented.str();
            stringFileSegmented.str("");

            reader->SetFileName(segmentedFile);
            reader->Update();

            InternalImageType::Pointer val = reader->GetOutput();
            double distance[4];
            distance[0] = 0.0;
            distance[1] = 0.0;
            distance[2] = 0.0;
            distance[3] = 0.0;
            distance[4] = 0.0;
            distance[5] = 0.0;
            distance[6] = 0.0;
            distance[7] = 0.0;
            distance[0] = calcDistance(val, x, y, 0);
            distance[1] = calcDistance(val, x, y, 1);
            distance[2] = calcDistance(val, x, y, 2);
            distance[3] = calcDistance(val, x, y, 3);
            distance[4] = calcDistance(val, x, y, 4);
            distance[5] = calcDistance(val, x, y, 5);
            distance[6] = calcDistance(val, x, y, 6);
            distance[7] = calcDistance(val, x, y, 7);
            extractValues<<"Pos 0 - Tickness"<<i+1<<": "<<distance[0]<<endl;
            extractValues<<"Pos 1 - Tickness"<<i+1<<": "<<distance[1]<<endl;
            extractValues<<"Pos 2 - Tickness"<<i+1<<": "<<distance[2]<<endl;
            extractValues<<"Pos 3 - Tickness"<<i+1<<": "<<distance[3]<<endl;
            extractValues<<"Pos 4 - Tickness"<<i+1<<": "<<distance[4]<<endl;
            extractValues<<"Pos 5 - Tickness"<<i+1<<": "<<distance[5]<<endl;
            extractValues<<"Pos 6 - Tickness"<<i+1<<": "<<distance[6]<<endl;
            extractValues<<"Pos 7 - Tickness"<<i+1<<": "<<distance[7]<<endl;
            extractValues<<"--------------------------------"<<endl;
        }
    }
    extractValues.close();
    endocardiumFile.close();

}

/**
 * @brief Distance::calcDistance
 * @param image
 * @param x
 * @param y
 * @param pos
 * @return double
 */
double Distance::calcDistance(ImageType2D::Pointer image, int x, int y, int pos){

    typedef unsigned int PixelType;

    typedef itk::Image< PixelType, 2> ImageType;

    bool isElement = false;
    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();

    itk::Point<int,2> p0;
    itk::Point<int,2> p1;
    if(pos == 0){
        for(int b= y; b > 0; b = b - 1){
            const ImageType::IndexType index = {{x,b}};
            if(((x) >= 0) && ((x) <= region[0]) && ((b) >= 0) && ((b) <= region[1])){
                if(image->GetPixel(index) > 0 && !isElement){
                    p0[0] = x;
                    p0[1] = b;
                    isElement = true;
                }
                if(isElement && image->GetPixel(index) > 0){
                    p1[0] = x;
                    p1[1] = b;
                }
            }
        }
    }

    if(pos == 1){
        for(int b= y; b<region[1]; b = b + 1){
            const ImageType::IndexType index = {{x,b}};
            if(((x) >= 0) && ((x) <= region[0]) && ((b) >= 0) && ((b) <= region[1])){
                if(image->GetPixel(index) > 0 && !isElement){
                    p0[0] = x;
                    p0[1] = b;
                    isElement = true;
                }
                if(isElement && image->GetPixel(index) > 0){
                    p1[0] = x;
                    p1[1] = b;
                }
            }
        }
    }

    if(pos == 2){
        for(int b= x; b > 0; b = b - 1){
            const ImageType::IndexType index = {{b,y}};
            if(((y) >= 0) && ((y) <= region[1]) && ((b) >= 0) && ((b) <= region[0])){
                if(image->GetPixel(index) > 0 && !isElement){
                    p0[0] = x;
                    p0[1] = b;
                    isElement = true;
                }
                if(isElement && image->GetPixel(index) > 0){
                    p1[0] = x;
                    p1[1] = b;
                }
            }
        }
    }

    if(pos == 3){
        for(int b= x; b<region[0]; b = b + 1){
            const ImageType::IndexType index = {{b,y}};
            if(((y) >= 0) && ((y) <= region[1]) && ((b) >= 0) && ((b) <= region[0])){
                if(image->GetPixel(index) > 0 && !isElement){
                    p0[0] = x;
                    p0[1] = b;
                    isElement = true;
                }
                if(isElement && image->GetPixel(index) > 0){
                    p1[0] = x;
                    p1[1] = b;
                }
            }
        }
    }

    if(pos == 4){
        int b = y;
        for(int a= x; a<region[0]; a = a + 1){
            if(b+1 < region[1]){
                const ImageType::IndexType index = {{a+1,b+1}};
                if(((a+1) >= 0) && ((a+1) <= region[0]) && ((b+1) >= 0) && ((b+1) <= region[1])){
                    if(image->GetPixel(index) > 0 && !isElement){
                        p0[0] = a;
                        p0[1] = b;
                        isElement = true;
                    }
                    if(isElement && image->GetPixel(index) > 0){
                        p1[0] = a;
                        p1[1] = b;
                    }
                }
            }
            b++;
        }
    }

    if(pos == 5){
        int b = y;
        for(int a= x; a>0; a = a - 1){
            if(b+1 < region[1]){
                const ImageType::IndexType index = {{a-1,b+1}};
                if(((a-1) >= 0) && ((a-1) <= region[0]) && ((b+1) >= 0) && ((b+1) <= region[1])){
                    if(image->GetPixel(index) > 0 && !isElement){
                        p0[0] = a;
                        p0[1] = b;
                        isElement = true;
                    }
                    if(isElement && image->GetPixel(index) > 0){
                        p1[0] = a;
                        p1[1] = b;
                    }
                }
            }
            b++;
        }
    }

    if(pos == 6){
        int a = x;
        for(int b=y; b>0; b = b - 1){
            if(a-1 > 0){
                const ImageType::IndexType index = {{a-1,b-1}};
                if(((a-1) >= 0) && ((a-1) <= region[0]) && ((b-1) >= 0) && ((b-1) <= region[1])){
                    if(image->GetPixel(index) > 0 && !isElement){
                        p0[0] = a;
                        p0[1] = b;
                        isElement = true;
                    }
                    if(isElement && image->GetPixel(index) > 0){
                        p1[0] = a;
                        p1[1] = b;
                    }
                }
            }
            a--;
        }
    }

    if(pos == 7){
        int a = x;
        for(int b=y; b>0; b = b - 1){
            if(a+1 <region[0]){
                const ImageType::IndexType index = {{a+1,b-1}};
                if(((a+1) >= 0) && ((a+1) <= region[0]) && ((b-1) >= 0) && ((b-1) <= region[1])){
                    if(image->GetPixel(index) > 0 && !isElement){
                        p0[0] = a;
                        p0[1] = b;
                        isElement = true;
                    }
                    if(isElement && image->GetPixel(index) > 0){
                        p1[0] = a;
                        p1[1] = b;
                    }
                }
            }
            a++;
        }
    }

    if(p0.EuclideanDistanceTo(p1) > 99999){
        return 0.0;
    }
    else{
        return p0.EuclideanDistanceTo(p1);
    }

}


