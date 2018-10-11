#include "kmeans.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkScalarImageKmeansImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionIterator.h"
#include <pwd.h>
#include <iostream>

using namespace std;

kmeans::kmeans()
{
    this->pw = getpwuid(getuid());

    this->homedir = this->pw->pw_dir;

    this->slices = "/temp/slices.txt";
    this->segmentedFinal = "/temp/segmentedFinal_";
    this->extractValues = "/temp/extractedValuesInternal.txt";
    this->pathKmeans = "/temp/classification/kmeans_";
    this->pathSegmentedFinal = this->homedir + this->segmentedFinal;
    this->pathExtractValues = this->homedir + this->extractValues;
    this->km = homedir + this->pathKmeans;

}

void kmeans::Execute(int first, int last){

    for(int i =first; i<last;i++){
        ImageType::Pointer image = ImageType::New();

        typedef itk::ScalarImageKmeansImageFilter< ImageType > KMeansFilterType;
        typedef itk::ImageFileReader<ImageType>  ReaderType;
        typedef itk::ImageFileWriter<ImageType> WriterType;

        KMeansFilterType::Pointer kmeansFilter = KMeansFilterType::New();
        string typeTiff = ".tif";
        stringstream segment;
        if(i<9)
            segment<<pathSegmentedFinal.c_str()<<"00"<<(i+1)<<typeTiff;
        if(i>=9 && i<99)
            segment<<pathSegmentedFinal.c_str()<<"0"<<(i+1)<<typeTiff;
        if(i>=99)
            segment<<pathSegmentedFinal.c_str()<<(i+1)<<typeTiff;
        string filenameSegmented = segment.str();
        segment.str("");

        typename ReaderType::Pointer readerFixed = ReaderType::New();
        readerFixed->SetFileName(filenameSegmented);
        readerFixed->Update();

        kmeansFilter->SetInput(readerFixed->GetOutput());
        kmeansFilter->SetUseNonContiguousLabels(true);
        kmeansFilter->AddClassWithInitialMean(0);
        kmeansFilter->AddClassWithInitialMean(64);
        kmeansFilter->AddClassWithInitialMean(128);
        kmeansFilter->AddClassWithInitialMean(192);
        kmeansFilter->Update();

        KMeansFilterType::ParametersType estimatedMeans = kmeansFilter->GetFinalMeans();

        const unsigned int numberOfClasses = estimatedMeans.Size();

        for(unsigned int i = 0 ; i < numberOfClasses ; ++i){
            std::cout << "cluster[" << i << "] ";
            std::cout << "    estimated mean : " << estimatedMeans[i] << std::endl;
        }

        typedef KMeansFilterType::OutputImageType  OutputImageType;

        typedef itk::RelabelComponentImageFilter<
                OutputImageType,
                OutputImageType > RelabelFilterType;

        RelabelFilterType::Pointer relabeler = RelabelFilterType::New();

        relabeler->SetInput( kmeansFilter->GetOutput() );

        typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleFilterType;
        RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
        rescaleFilter->SetInput(relabeler->GetOutput());
        rescaleFilter->SetOutputMinimum(0);
        rescaleFilter->SetOutputMaximum(255);

        typedef RelabelFilterType::ObjectSizeInPixelsContainerType SizesType;

        const SizesType &  sizes = relabeler->GetSizeOfObjectsInPixels();

        SizesType::const_iterator sizeItr = sizes.begin();
        SizesType::const_iterator sizeEnd = sizes.end();

        std::cout << "Number of pixels per class " << std::endl;
        unsigned int kclass = 0;
        while( sizeItr != sizeEnd ){
            std::cout << "Class " << kclass << " = " << *sizeItr << std::endl;
            ++kclass;
            ++sizeItr;
        }

        typename WriterType::Pointer writer = WriterType::New();

        stringstream stringFile;
        stringFile<<km<<(i+1)<<typeTiff;

        string File = stringFile.str();
        stringFile.str("");

        writer->SetFileName(File);
        writer->SetInput( rescaleFilter->GetOutput() );
        writer->Update();
    }

}
