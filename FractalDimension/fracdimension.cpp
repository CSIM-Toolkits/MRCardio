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

using namespace std;

fracdimension::fracdimension()
{
}

//double fracdimension::GetDimension(ImageType::Pointer Image){

//    double dim = 0.0;
//    int r = 2;
//    int cont_box[6];
//    for(int i = 0;i<5;i++)
//        cont_box[i]=0;
//    double la[5];
//    double sum_la = 0.0;
//    typedef unsigned char PixelType;
//    typedef itk::Image<PixelType,3> ImageType;
//    typedef itk::ShapedNeighborhoodIterator<ImageType> IteratorType;
//    typedef itk::ConstNeighborhoodIterator<ImageType> ShapedNeighborhoodIteratorType;
//    typedef itk::FlatStructuringElement<3> FlatStructuringElementType;
//    FlatStructuringElementType::RadiusType r1;
//    r1.Fill(r);
//    //ShapedNeighborhoodIteratorType::RadiusType radius;
//    //radius.Fill(3);
//    itk::Size<3> radius;
//    radius.Fill(2);

//    IteratorType iterator(radius,Image,Image->GetLargestPossibleRegion());
//    /*
//    IteratorType::OffsetType top = {{0,-1}};
//    iterator.ActivateOffset((top));
//    IteratorType::OffsetType bottom = {{0,1}};
//    iterator.ActivateOffset(bottom);
//    IteratorType::OffsetType left = {{-1,0}};
//    iterator.ActivateOffset(left);
//    IteratorType::OffsetType right = {{1,0}};
//    iterator.ActivateOffset((right));

//    IteratorType::OffsetType topr = {{1,-1}};
//    iterator.ActivateOffset((topr));
//    IteratorType::OffsetType bottomr = {{1,1}};
//    iterator.ActivateOffset(bottomr);
//    IteratorType::OffsetType topl = {{-1,-1}};
//    iterator.ActivateOffset(topl);
//    IteratorType::OffsetType bottoml = {{-1,1}};
//    iterator.ActivateOffset((bottoml));
//    IteratorType::OffsetType center = {{0,0}};
//    iterator.ActivateOffset((center));
//    */

//    bool flag = false;
//    IteratorType::OffsetType off;
//    for(int i = 0; i<3; i++){
//        for (int y = -r; y <= r; y++){
//            for (int x = -r; x <= r; x++){
//                //for(int z = -r; z <= r; z++){

//                off[0] = x;
//                off[1] = y;
//                //off[2] = 1;
//                iterator.ActivateOffset(off);
//                //}
//            }
//        }
//        //IteratorType::ConstIterator es = iterator.Begin();
//        for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator){
//            IteratorType::ConstIterator es = iterator.Begin();
//            flag = false;
//            while(!es.IsAtEnd()){
//                if(es.Get()!=0){
//                    flag = true;
//                }
//                es++;
//            }
//            if(flag==true)
//                cont_box[i]++;
//        }
//        r = r * 2;
//    }
//    r = 2;
//    for(int k=0;k<3;k++){
//        la[k] = log(cont_box[k])/log(r);
//        cout<<"dim "<<k<<" = "<<la[k]<<endl;
//        sum_la = sum_la + la[k];
//        r = r * 2;
//    }
//    dim = -(sum_la)/3;

//    return dim;
//}

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

    typedef itk::BinaryBallStructuringElement< PixelType, 3>
            StructuringElementType;
    StructuringElementType structuringElement;
    structuringElement.SetRadius(1);
    structuringElement.CreateStructuringElement();

    IteratorType iterator(structuringElement.GetRadius(), image, image->GetLargestPossibleRegion());
    iterator.CreateActiveListFromNeighborhood(structuringElement);

    const ImageType::SizeType region = image->GetLargestPossibleRegion().GetSize();
    double M = region[0];
    double r = 1/M;

    IteratorType::IndexListType indexList = iterator.GetActiveIndexList();
    IteratorType::IndexListType::const_iterator
            listIterator = indexList.begin();

    // Note that ZeroFluxNeumannBoundaryCondition is used by default so even
    // pixels outside of the image will have valid values (equivalent to
    // their neighbors just inside the image)
    for( iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator )
    {
        IteratorType::ConstIterator ci = iterator.Begin();

        while( !ci.IsAtEnd() )
        {
            isZero = false;
            isElement = false;
            if(ci.Get() == 0){
                isZero = true;
            }
            if(ci.Get() !=0){
                isElement = true;
            }
            std::cout << "Centered at " << iterator.GetIndex() << std::endl;
            std::cout << "Neighborhood index " << ci.GetNeighborhoodIndex()
              << " is offset " << ci.GetNeighborhoodOffset()
              << " and has value " << ci.Get()
              << " The real index is "
              << iterator.GetIndex() + ci.GetNeighborhoodOffset()
              << std::endl;
            ++ci;
        }
        if(isZero == 1 && isElement == 1){
            cont_A++;
        }
    }
    cout<<"NR: "<<cont_A<<endl;
    cout<<"LOG NR: "<<log(cont_A)<<endl;
    cout<<"r: "<<r<<endl;
    cout<<"1/r: "<<1/r<<endl;
    cout<<"log(1/r): "<<log(1/r)<<endl;
    dim = ((log(cont_A)) / (log(1/r)));
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

void fracdimension::GradientMagnitude(ImageType::Pointer Image, const char* volume_out){
    typedef itk::GradientMagnitudeImageFilter<
            ImageType, ImageType >  filterType;

    // Create and setup a gradient filter
    filterType::Pointer gradientFilter = filterType::New();
    gradientFilter->SetInput(Image);
    gradientFilter->Update();
    typedef itk::BinaryThresholdImageFilter<
            ImageType, ImageType>  FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(gradientFilter->GetOutput());
    //filter->SetUseImageSpacing( useImageSpacing );
    filter->SetLowerThreshold(10);
    filter->SetUpperThreshold(255);
    filter->SetOutsideValue(0);
    filter->SetInsideValue(255);
    filter->Update();

    edge = filter->GetOutput();

    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(volume_out);
    writer->SetInput(filter->GetOutput());
    writer->Update();

}

void fracdimension::MorphologicalGradient(ImageType::Pointer Image){
    ImageType::Pointer edge;
    unsigned int radius = 3;
    typedef itk::FlatStructuringElement<3> FlatStructuringElementType;
    FlatStructuringElementType::RadiusType flat;
    flat.Fill(radius);

    FlatStructuringElementType structuringElement = FlatStructuringElementType::Box(flat);

    typedef itk::BinaryErodeImageFilter <ImageType, ImageType, FlatStructuringElementType>
            BinaryErodeImageFilterType;

    BinaryErodeImageFilterType::Pointer erodeFilter
            = BinaryErodeImageFilterType::New();
    erodeFilter->SetInput(Image);
    erodeFilter->SetKernel(structuringElement);
    //erodeFilter->SetErodeValue(255);

    typedef itk::BinaryDilateImageFilter <ImageType, ImageType, FlatStructuringElementType>
            BinaryDilateImageFilterType;

    BinaryDilateImageFilterType::Pointer dilateFilter
            = BinaryDilateImageFilterType::New();
    dilateFilter->SetInput(Image);
    dilateFilter->SetKernel(structuringElement);
    dilateFilter->SetDilateValue(255);

    typedef itk::SubtractImageFilter<ImageType,ImageType,ImageType> SubtractImageFilterType;
    SubtractImageFilterType::Pointer subtract = SubtractImageFilterType::New();
    subtract->SetInput1(dilateFilter->GetOutput());
    subtract->SetInput2(erodeFilter->GetOutput());
    //subtract->Update();

    //edge = subtract->GetOutput();

    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName("/home/gustavo/temp/outputGradient.tif");
    writer->SetInput(subtract->GetOutput());
    writer->Update();

}
