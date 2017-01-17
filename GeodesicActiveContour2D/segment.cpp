#include "segment.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include <iostream>
#include "QString"
#include "itkOrientImageFilter.h"
#include <fstream>
#include <string>
#include <sstream>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkAndImageFilter.h>
#include "utils.h"

using namespace std;

segment::segment()
{
}

void segment::InternalEC(int first,int last, double sigma, double sig_min, double sig_max, double propagation, double curvature, double advection, double rms, int iterations, double timestep, int it_dif, double conductance, double alpha, double beta, double distance){
    ifstream myfile ("/home/gustavo/temp/endocardium.txt");
    for(int i =first; i<last;i++){
        typedef   float           InternalPixelType;
        const     unsigned int    Dimension = 2;

        typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;

        typedef unsigned char                            OutputPixelType;
        typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
        typedef itk::BinaryThresholdImageFilter<
                InternalImageType,
                OutputImageType    >       ThresholdingFilterType;
        ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
        thresholder->SetLowerThreshold( -1000.0 );
        thresholder->SetUpperThreshold(     0.0 );
        thresholder->SetOutsideValue(  0  );
        thresholder->SetInsideValue(  255 );

        typedef  itk::ImageFileReader< InternalImageType > ReaderType;
        typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;
        ReaderType::Pointer reader = ReaderType::New();
        WriterType::Pointer writer = WriterType::New();
        //WriterType::Pointer writer_out = WriterType::New();
        stringstream ss;

        string name = "/home/gustavo/temp/cine_";
        string type = ".tif";

        ss<<name<<(i+1)<<type;

        string filename = ss.str();
        ss.str("");

        reader->SetFileName(filename);
        reader->Update();

        InternalImageType::Pointer val = reader->GetOutput();

        typedef itk::RescaleIntensityImageFilter<
                InternalImageType,
                OutputImageType >   CastFilterType;

        typedef   itk::CurvatureAnisotropicDiffusionImageFilter<
                InternalImageType,
                InternalImageType >  SmoothingFilterType;
        SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();

        typedef   itk::GradientMagnitudeRecursiveGaussianImageFilter<
                InternalImageType,
                InternalImageType >  GradientFilterType;
        typedef   itk::SigmoidImageFilter<
                InternalImageType,
                InternalImageType >  SigmoidFilterType;
        GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
        SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();

        sigmoid->SetOutputMinimum(sig_min);
        sigmoid->SetOutputMaximum(sig_max);

        typedef  itk::FastMarchingImageFilter<
                InternalImageType,
                InternalImageType >    FastMarchingFilterType;

        FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
        //const InternalImageType * inputImage = reader->GetOutput();
        const InternalImageType * inputImage = val;
        fastMarching->SetOutputRegion( inputImage->GetBufferedRegion() );
        fastMarching->SetOutputSpacing( inputImage->GetSpacing() );
        fastMarching->SetOutputOrigin( inputImage->GetOrigin() );
        fastMarching->SetOutputDirection( inputImage->GetDirection() );

        typedef  itk::GeodesicActiveContourLevelSetImageFilter< InternalImageType,
                InternalImageType >    GeodesicActiveContourFilterType;
        GeodesicActiveContourFilterType::Pointer geodesicActiveContour =
                GeodesicActiveContourFilterType::New();

        const double propagationScaling = propagation;
        //  Software Guide : BeginCodeSnippet
        geodesicActiveContour->SetPropagationScaling( propagationScaling );
        geodesicActiveContour->SetCurvatureScaling(curvature);
        geodesicActiveContour->SetAdvectionScaling(advection);

        geodesicActiveContour->SetMaximumRMSError(rms);
        geodesicActiveContour->SetNumberOfIterations(iterations);

        smoothing->SetInput(val);
        gradientMagnitude->SetInput( smoothing->GetOutput() );
        sigmoid->SetInput( gradientMagnitude->GetOutput() );
        geodesicActiveContour->SetInput(  fastMarching->GetOutput() );
        geodesicActiveContour->SetFeatureImage( sigmoid->GetOutput() );
        thresholder->SetInput( geodesicActiveContour->GetOutput() );
        //ImageOut = thresholder->GetOutput();
        //writer->SetInput( thresholder->GetOutput() );

        smoothing->SetTimeStep(timestep);
        smoothing->SetNumberOfIterations(it_dif);
        smoothing->SetConductanceParameter(conductance);

        const double sig = sigma;
        gradientMagnitude->SetSigma(sig);

        sigmoid->SetAlpha( alpha );
        sigmoid->SetBeta(  beta  );

        string line;
        int x;
        int y;

        if (myfile.is_open())
        {
            getline (myfile,line);
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

        typedef FastMarchingFilterType::NodeContainer  NodeContainer;
        typedef FastMarchingFilterType::NodeType       NodeType;
        NodeContainer::Pointer seeds = NodeContainer::New();
        InternalImageType::IndexType  seedPosition;
        seedPosition.SetElement(0,x);
        seedPosition.SetElement(1,y);
        //seedPosition.SetElement(2,(int)z);
        cout<<"X, Y = "<<x<<" "<<y<<endl;
        NodeType node;
        const double seedValue = - distance;
        node.SetValue( seedValue );
        node.SetIndex( seedPosition );

        seeds->Initialize();
        seeds->InsertElement( 0, node );

        fastMarching->SetTrialPoints(  seeds  );

        fastMarching->SetSpeedConstant( 1.0 );

        CastFilterType::Pointer caster1 = CastFilterType::New();
        CastFilterType::Pointer caster2 = CastFilterType::New();
        CastFilterType::Pointer caster3 = CastFilterType::New();
        CastFilterType::Pointer caster4 = CastFilterType::New();
        WriterType::Pointer writer1 = WriterType::New();
        WriterType::Pointer writer2 = WriterType::New();
        WriterType::Pointer writer3 = WriterType::New();
        WriterType::Pointer writer4 = WriterType::New();
        caster1->SetInput( smoothing->GetOutput() );
        writer1->SetInput( caster1->GetOutput() );
        writer1->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput1.tif");
        caster1->SetOutputMinimum(   0 );
        caster1->SetOutputMaximum( 255 );
        writer1->Update();
        caster2->SetInput( gradientMagnitude->GetOutput() );
        writer2->SetInput( caster2->GetOutput() );
        writer2->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput2.tif");
        caster2->SetOutputMinimum(   0 );
        caster2->SetOutputMaximum( 255 );
        writer2->Update();
        caster3->SetInput( sigmoid->GetOutput() );
        writer3->SetInput( caster3->GetOutput() );
        writer3->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput3.tif");
        caster3->SetOutputMinimum(   0 );
        caster3->SetOutputMaximum( 255 );
        writer3->Update();
        caster4->SetInput( fastMarching->GetOutput() );
        writer4->SetInput( caster4->GetOutput() );
        writer4->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput4.tif");
        caster4->SetOutputMinimum(   0 );
        caster4->SetOutputMaximum( 255 );

        fastMarching->SetOutputSize(
                    reader->GetOutput()->GetBufferedRegion().GetSize() );
        reader->Update();

        stringstream ss2;

        string name2 = "/home/gustavo/temp/segmented_";
        string type2 = ".tif";

        if(i<9)
            ss2<<name2<<"00"<<(i+1)<<type2;
        if(i>=9 && i<99)
            ss2<<name2<<"0"<<(i+1)<<type2;
        if(i>=99)
            ss2<<name2<<(i+1)<<type2;
        string filename2 = ss2.str();
        ss2.str("");

        typedef itk::GradientMagnitudeImageFilter<OutputImageType, OutputImageType >  GradientType;
        GradientType::Pointer gradient = GradientType::New();
        gradient->SetInput(thresholder->GetOutput());
        gradient->Update();

        try
        {
            writer->SetFileName(filename2);
            writer->SetInput( gradient->GetOutput() );
            writer->Update();
        }
        catch( itk::ExceptionObject & excep )
        {
            std::cerr << "Exception caught !" << std::endl;
            std::cerr << excep << std::endl;
        }

        std::cout << std::endl;
        std::cout << "Max. no. iterations: " << geodesicActiveContour->GetNumberOfIterations() << std::endl;
        std::cout << "Max. RMS error: " << geodesicActiveContour->GetMaximumRMSError() << std::endl;
        std::cout << std::endl;
        std::cout << "No. elpased iterations: " << geodesicActiveContour->GetElapsedIterations() << std::endl;
        std::cout << "RMS change: " << geodesicActiveContour->GetRMSChange() << std::endl;
        writer4->Update();

        typedef itk::ImageFileWriter< InternalImageType > InternalWriterType;
        InternalWriterType::Pointer mapWriter = InternalWriterType::New();
        mapWriter->SetInput( fastMarching->GetOutput() );
        mapWriter->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput4.mha");
        mapWriter->Update();
        InternalWriterType::Pointer speedWriter = InternalWriterType::New();
        speedWriter->SetInput( sigmoid->GetOutput() );
        speedWriter->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput3.mha");
        speedWriter->Update();
        InternalWriterType::Pointer gradientWriter = InternalWriterType::New();
        gradientWriter->SetInput( gradientMagnitude->GetOutput() );
        gradientWriter->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput2.mha");
        gradientWriter->Update();
    }
    myfile.close();

}

void segment::MyocardiumEC(int first,int last, double sigma, double sig_min, double sig_max, double propagation, double curvature, double advection, double rms, int iterations, double timestep, int it_dif, double conductance, double alpha, double beta, double distance){
    ifstream myfile ("/home/gustavo/temp/endocardium.txt");
    ifstream myfile2 ("/home/gustavo/temp/radius.txt");
    for(int i =first; i<last;i++){
        typedef   float           InternalPixelType;
        const     unsigned int    Dimension = 2;

        typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;

        typedef unsigned char                            OutputPixelType;
        typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
        typedef itk::BinaryThresholdImageFilter<
                InternalImageType,
                OutputImageType    >       ThresholdingFilterType;
        ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
        thresholder->SetLowerThreshold( -1000.0 );
        thresholder->SetUpperThreshold(     0.0 );
        thresholder->SetOutsideValue(  0  );
        thresholder->SetInsideValue(  255 );

        typedef  itk::ImageFileReader< InternalImageType > ReaderType;
        typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;
        ReaderType::Pointer reader = ReaderType::New();
        WriterType::Pointer writer = WriterType::New();
        //WriterType::Pointer writer_out = WriterType::New();
        stringstream ss;

        string name = "/home/gustavo/temp/cine_";
        string type = ".tif";

        ss<<name<<(i+1)<<type;

        string filename = ss.str();
        ss.str("");

        reader->SetFileName(filename);
        reader->Update();

        InternalImageType::Pointer val = reader->GetOutput();

        typedef itk::RescaleIntensityImageFilter<
                InternalImageType,
                OutputImageType >   CastFilterType;

        typedef   itk::CurvatureAnisotropicDiffusionImageFilter<
                InternalImageType,
                InternalImageType >  SmoothingFilterType;
        SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();

        typedef   itk::GradientMagnitudeRecursiveGaussianImageFilter<
                InternalImageType,
                InternalImageType >  GradientFilterType;
        typedef   itk::SigmoidImageFilter<
                InternalImageType,
                InternalImageType >  SigmoidFilterType;
        GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
        SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();

        sigmoid->SetOutputMinimum(sig_min);
        sigmoid->SetOutputMaximum(sig_max);

        typedef  itk::FastMarchingImageFilter<
                InternalImageType,
                InternalImageType >    FastMarchingFilterType;

        FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
        //const InternalImageType * inputImage = reader->GetOutput();
        const InternalImageType * inputImage = val;
        fastMarching->SetOutputRegion( inputImage->GetBufferedRegion() );
        fastMarching->SetOutputSpacing( inputImage->GetSpacing() );
        fastMarching->SetOutputOrigin( inputImage->GetOrigin() );
        fastMarching->SetOutputDirection( inputImage->GetDirection() );

        typedef  itk::GeodesicActiveContourLevelSetImageFilter< InternalImageType,
                InternalImageType >    GeodesicActiveContourFilterType;
        GeodesicActiveContourFilterType::Pointer geodesicActiveContour =
                GeodesicActiveContourFilterType::New();

        const double propagationScaling = propagation;
        //  Software Guide : BeginCodeSnippet
        geodesicActiveContour->SetPropagationScaling( propagationScaling );
        geodesicActiveContour->SetCurvatureScaling(curvature);
        geodesicActiveContour->SetAdvectionScaling(advection);

        geodesicActiveContour->SetMaximumRMSError(rms);
        geodesicActiveContour->SetNumberOfIterations(iterations);

        smoothing->SetInput(val);
        gradientMagnitude->SetInput( smoothing->GetOutput() );
        sigmoid->SetInput( gradientMagnitude->GetOutput() );
        geodesicActiveContour->SetInput(  fastMarching->GetOutput() );
        geodesicActiveContour->SetFeatureImage( sigmoid->GetOutput() );
        thresholder->SetInput( geodesicActiveContour->GetOutput() );
        //ImageOut = thresholder->GetOutput();
        //writer->SetInput( thresholder->GetOutput() );

        smoothing->SetTimeStep(timestep);
        smoothing->SetNumberOfIterations(it_dif);
        smoothing->SetConductanceParameter(conductance);

        const double sig = sigma;
        gradientMagnitude->SetSigma(sig);

        sigmoid->SetAlpha( alpha );
        sigmoid->SetBeta(  beta  );

        string line;
        string line2;
        int x;
        int y;
        int radius;
        if (myfile.is_open())
        {
            getline (myfile,line);
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
        if (myfile2.is_open())
        {
            getline(myfile2,line2);
            radius = atoi(line2.c_str());
        }

        typedef FastMarchingFilterType::NodeContainer  NodeContainer;
        typedef FastMarchingFilterType::NodeType       NodeType;
        NodeContainer::Pointer seeds = NodeContainer::New();
        InternalImageType::IndexType  seedPosition;
        seedPosition.SetElement(0,x);
        seedPosition.SetElement(1,(y+(radius-2)));
        //seedPosition.SetElement(2,(int)z);
        cout<<"X, Y = "<<x<<" "<<y<<endl;
        NodeType node;
        const double seedValue = - distance;
        node.SetValue( seedValue );
        node.SetIndex( seedPosition );

        seeds->Initialize();
        seeds->InsertElement( 0, node );

        fastMarching->SetTrialPoints(  seeds  );

        fastMarching->SetSpeedConstant( 1.0 );

        CastFilterType::Pointer caster1 = CastFilterType::New();
        CastFilterType::Pointer caster2 = CastFilterType::New();
        CastFilterType::Pointer caster3 = CastFilterType::New();
        CastFilterType::Pointer caster4 = CastFilterType::New();
        WriterType::Pointer writer1 = WriterType::New();
        WriterType::Pointer writer2 = WriterType::New();
        WriterType::Pointer writer3 = WriterType::New();
        WriterType::Pointer writer4 = WriterType::New();
        caster1->SetInput( smoothing->GetOutput() );
        writer1->SetInput( caster1->GetOutput() );
        writer1->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput1_2.tif");
        caster1->SetOutputMinimum(   0 );
        caster1->SetOutputMaximum( 255 );
        writer1->Update();
        caster2->SetInput( gradientMagnitude->GetOutput() );
        writer2->SetInput( caster2->GetOutput() );
        writer2->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput2_2.tif");
        caster2->SetOutputMinimum(   0 );
        caster2->SetOutputMaximum( 255 );
        writer2->Update();
        caster3->SetInput( sigmoid->GetOutput() );
        writer3->SetInput( caster3->GetOutput() );
        writer3->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput3_2.tif");
        caster3->SetOutputMinimum(   0 );
        caster3->SetOutputMaximum( 255 );
        writer3->Update();
        caster4->SetInput( fastMarching->GetOutput() );
        writer4->SetInput( caster4->GetOutput() );
        writer4->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput4_2.tif");
        caster4->SetOutputMinimum(   0 );
        caster4->SetOutputMaximum( 255 );

        fastMarching->SetOutputSize(
                    reader->GetOutput()->GetBufferedRegion().GetSize() );
        reader->Update();

        std::cout << std::endl;
        std::cout << "Max. no. iterations: " << geodesicActiveContour->GetNumberOfIterations() << std::endl;
        std::cout << "Max. RMS error: " << geodesicActiveContour->GetMaximumRMSError() << std::endl;
        std::cout << std::endl;
        std::cout << "No. elpased iterations: " << geodesicActiveContour->GetElapsedIterations() << std::endl;
        std::cout << "RMS change: " << geodesicActiveContour->GetRMSChange() << std::endl;
        writer4->Update();

        typedef itk::ImageFileWriter< InternalImageType > InternalWriterType;
        InternalWriterType::Pointer mapWriter = InternalWriterType::New();
        mapWriter->SetInput( fastMarching->GetOutput() );
        mapWriter->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput4_2.mha");
        mapWriter->Update();
        InternalWriterType::Pointer speedWriter = InternalWriterType::New();
        speedWriter->SetInput( sigmoid->GetOutput() );
        speedWriter->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput3_2.mha");
        speedWriter->Update();
        InternalWriterType::Pointer gradientWriter = InternalWriterType::New();
        gradientWriter->SetInput( gradientMagnitude->GetOutput() );
        gradientWriter->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput2_2.mha");
        gradientWriter->Update();

        CastFilterType::Pointer caster5 = CastFilterType::New();
        caster5->SetInput( reader->GetOutput() );
        caster5->SetOutputMinimum(   0 );
        caster5->SetOutputMaximum( 255 );
        caster5->Update();
        typedef itk::AndImageFilter <OutputImageType>  AndImageFilterType;
        AndImageFilterType::Pointer andFilter
                = AndImageFilterType::New();
        andFilter->SetInput(0, caster5->GetOutput());
        andFilter->SetInput(1, thresholder->GetOutput());
        andFilter->Update();

        stringstream ss2;

        string name2 = "/home/gustavo/temp/segmented2_";
        string type2 = ".tif";

        if(i<9)
            ss2<<name2<<"00"<<(i+1)<<type2;
        if(i>=9 && i<99)
            void MyocardiumEC(int first,int last,double sigma,double sig_min, double sig_max, double propagation, double curvature, double advection, double rms, int iterations, double timestep, int it_dif, double conductance, double alpha, double meta, double distance);
            ss2<<name2<<"0"<<(i+1)<<type2;
        if(i>=99)
            ss2<<name2<<(i+1)<<type2;
        string filename2 = ss2.str();
        ss2.str("");

//        typedef itk::GradientMagnitudeImageFilter<OutputImageType, OutputImageType >  GradientType;
//        GradientType::Pointer gradient = GradientType::New();
//        gradient->SetInput(thresholder->GetOutput());
//        gradient->Update();

        try
        {
            writer->SetFileName(filename2);
            writer->SetInput(andFilter->GetOutput() );
            writer->Update();
        }
        catch( itk::ExceptionObject & excep )
        {
            std::cerr << "Exception caught !" << std::endl;
            std::cerr << excep << std::endl;
        }

    }
    myfile.close();

}

void segment::InternalELV(int first,int last, double sigma, double sig_min, double sig_max, double propagation, double curvature, double advection, double rms, int iterations, double timestep, int it_dif, double conductance, double alpha, double beta, double distance){
    ifstream myfile ("/home/gustavo/temp/endocardium.txt");
    for(int i =first; i<last;i++){
        typedef   float           InternalPixelType;
        const     unsigned int    Dimension = 2;

        typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;

        typedef unsigned char                            OutputPixelType;
        typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
        typedef itk::BinaryThresholdImageFilter<
                InternalImageType,
                OutputImageType    >       ThresholdingFilterType;
        ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
        thresholder->SetLowerThreshold( -1000.0 );
        thresholder->SetUpperThreshold(     0.0 );
        thresholder->SetOutsideValue(  0  );
        thresholder->SetInsideValue(  255 );

        typedef  itk::ImageFileReader< InternalImageType > ReaderType;
        typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;
        ReaderType::Pointer reader = ReaderType::New();
        WriterType::Pointer writer = WriterType::New();
        //WriterType::Pointer writer_out = WriterType::New();
        stringstream ss;

        string name = "/home/gustavo/temp/cine_";
        string type = ".tif";

        ss<<name<<(i+1)<<type;

        string filename = ss.str();
        ss.str("");

        reader->SetFileName(filename);
        reader->Update();

        InternalImageType::Pointer val = reader->GetOutput();

        typedef itk::RescaleIntensityImageFilter<
                InternalImageType,
                OutputImageType >   CastFilterType;

        typedef   itk::CurvatureAnisotropicDiffusionImageFilter<
                InternalImageType,
                InternalImageType >  SmoothingFilterType;
        SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();

        typedef   itk::GradientMagnitudeRecursiveGaussianImageFilter<
                InternalImageType,
                InternalImageType >  GradientFilterType;
        typedef   itk::SigmoidImageFilter<
                InternalImageType,
                InternalImageType >  SigmoidFilterType;
        GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
        SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();

        sigmoid->SetOutputMinimum(sig_min);
        sigmoid->SetOutputMaximum(sig_max);

        typedef  itk::FastMarchingImageFilter<
                InternalImageType,
                InternalImageType >    FastMarchingFilterType;

        FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
        //const InternalImageType * inputImage = reader->GetOutput();
        const InternalImageType * inputImage = val;
        fastMarching->SetOutputRegion( inputImage->GetBufferedRegion() );
        fastMarching->SetOutputSpacing( inputImage->GetSpacing() );
        fastMarching->SetOutputOrigin( inputImage->GetOrigin() );
        fastMarching->SetOutputDirection( inputImage->GetDirection() );

        typedef  itk::GeodesicActiveContourLevelSetImageFilter< InternalImageType,
                InternalImageType >    GeodesicActiveContourFilterType;
        GeodesicActiveContourFilterType::Pointer geodesicActiveContour =
                GeodesicActiveContourFilterType::New();

        const double propagationScaling = propagation;
        //  Software Guide : BeginCodeSnippet
        geodesicActiveContour->SetPropagationScaling( propagationScaling );
        geodesicActiveContour->SetCurvatureScaling(curvature);
        geodesicActiveContour->SetAdvectionScaling(advection);

        geodesicActiveContour->SetMaximumRMSError(rms);
        geodesicActiveContour->SetNumberOfIterations(iterations);

        smoothing->SetInput(val);
        gradientMagnitude->SetInput( smoothing->GetOutput() );
        sigmoid->SetInput( gradientMagnitude->GetOutput() );
        geodesicActiveContour->SetInput(  fastMarching->GetOutput() );
        geodesicActiveContour->SetFeatureImage( sigmoid->GetOutput() );
        thresholder->SetInput( geodesicActiveContour->GetOutput() );
        //ImageOut = thresholder->GetOutput();
        //writer->SetInput( thresholder->GetOutput() );

        smoothing->SetTimeStep(timestep);
        smoothing->SetNumberOfIterations(it_dif);
        smoothing->SetConductanceParameter(conductance);

        const double sig = sigma;
        gradientMagnitude->SetSigma(sig);

        sigmoid->SetAlpha( alpha );
        sigmoid->SetBeta(  beta  );

        string line;
        int x;
        int y;

        if (myfile.is_open())
        {
            getline (myfile,line);
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

        typedef FastMarchingFilterType::NodeContainer  NodeContainer;
        typedef FastMarchingFilterType::NodeType       NodeType;
        NodeContainer::Pointer seeds = NodeContainer::New();
        InternalImageType::IndexType  seedPosition;
        seedPosition.SetElement(0,x);
        seedPosition.SetElement(1,y);
        //seedPosition.SetElement(2,(int)z);
        cout<<"X, Y = "<<x<<" "<<y<<endl;
        NodeType node;
        const double seedValue = - distance;
        node.SetValue( seedValue );
        node.SetIndex( seedPosition );

        seeds->Initialize();
        seeds->InsertElement( 0, node );

        fastMarching->SetTrialPoints(  seeds  );

        fastMarching->SetSpeedConstant( 1.0 );

        CastFilterType::Pointer caster1 = CastFilterType::New();
        CastFilterType::Pointer caster2 = CastFilterType::New();
        CastFilterType::Pointer caster3 = CastFilterType::New();
        CastFilterType::Pointer caster4 = CastFilterType::New();
        WriterType::Pointer writer1 = WriterType::New();
        WriterType::Pointer writer2 = WriterType::New();
        WriterType::Pointer writer3 = WriterType::New();
        WriterType::Pointer writer4 = WriterType::New();
        caster1->SetInput( smoothing->GetOutput() );
        writer1->SetInput( caster1->GetOutput() );
        writer1->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput1.tif");
        caster1->SetOutputMinimum(   0 );
        caster1->SetOutputMaximum( 255 );
        writer1->Update();
        caster2->SetInput( gradientMagnitude->GetOutput() );
        writer2->SetInput( caster2->GetOutput() );
        writer2->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput2.tif");
        caster2->SetOutputMinimum(   0 );
        caster2->SetOutputMaximum( 255 );
        writer2->Update();
        caster3->SetInput( sigmoid->GetOutput() );
        writer3->SetInput( caster3->GetOutput() );
        writer3->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput3.tif");
        caster3->SetOutputMinimum(   0 );
        caster3->SetOutputMaximum( 255 );
        writer3->Update();
        caster4->SetInput( fastMarching->GetOutput() );
        writer4->SetInput( caster4->GetOutput() );
        writer4->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput4.tif");
        caster4->SetOutputMinimum(   0 );
        caster4->SetOutputMaximum( 255 );

        fastMarching->SetOutputSize(
                    reader->GetOutput()->GetBufferedRegion().GetSize() );
        reader->Update();

        stringstream ss2;

        string name2 = "/home/gustavo/temp/segmented_";
        string type2 = ".tif";

        if(i<9)
            ss2<<name2<<"00"<<(i+1)<<type2;
        if(i>=9 && i<99)
            ss2<<name2<<"0"<<(i+1)<<type2;
        if(i>=99)
            ss2<<name2<<(i+1)<<type2;
        string filename2 = ss2.str();
        ss2.str("");

        typedef itk::GradientMagnitudeImageFilter<OutputImageType, OutputImageType >  GradientType;
        GradientType::Pointer gradient = GradientType::New();
        gradient->SetInput(thresholder->GetOutput());
        gradient->Update();

        try
        {
            writer->SetFileName(filename2);
            writer->SetInput( thresholder->GetOutput() );
            writer->Update();
        }
        catch( itk::ExceptionObject & excep )
        {
            std::cerr << "Exception caught !" << std::endl;
            std::cerr << excep << std::endl;
        }

        std::cout << std::endl;
        std::cout << "Max. no. iterations: " << geodesicActiveContour->GetNumberOfIterations() << std::endl;
        std::cout << "Max. RMS error: " << geodesicActiveContour->GetMaximumRMSError() << std::endl;
        std::cout << std::endl;
        std::cout << "No. elpased iterations: " << geodesicActiveContour->GetElapsedIterations() << std::endl;
        std::cout << "RMS change: " << geodesicActiveContour->GetRMSChange() << std::endl;
        writer4->Update();

        typedef itk::ImageFileWriter< InternalImageType > InternalWriterType;
        InternalWriterType::Pointer mapWriter = InternalWriterType::New();
        mapWriter->SetInput( fastMarching->GetOutput() );
        mapWriter->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput4.mha");
        mapWriter->Update();
        InternalWriterType::Pointer speedWriter = InternalWriterType::New();
        speedWriter->SetInput( sigmoid->GetOutput() );
        speedWriter->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput3.mha");
        speedWriter->Update();
        InternalWriterType::Pointer gradientWriter = InternalWriterType::New();
        gradientWriter->SetInput( gradientMagnitude->GetOutput() );
        gradientWriter->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput2.mha");
        gradientWriter->Update();
    }
    myfile.close();

}

void segment::MyocardiumELV(int first,int last, double sigma, double sig_min, double sig_max, double propagation, double curvature, double advection, double rms, int iterations, double timestep, int it_dif, double conductance, double alpha, double beta, double distance){
    ifstream myfile ("/home/gustavo/temp/endocardium.txt");
    ifstream myfile2 ("/home/gustavo/temp/radius.txt");
    for(int i =first; i<last;i++){
        typedef   float           InternalPixelType;
        const     unsigned int    Dimension = 2;

        typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;

        typedef unsigned char                            OutputPixelType;
        typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
        typedef itk::BinaryThresholdImageFilter<
                InternalImageType,
                OutputImageType    >       ThresholdingFilterType;
        ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
        thresholder->SetLowerThreshold( -1000.0 );
        thresholder->SetUpperThreshold(     0.0 );
        thresholder->SetOutsideValue(  0  );
        thresholder->SetInsideValue(  255 );

        typedef  itk::ImageFileReader< InternalImageType > ReaderType;
        typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;
        ReaderType::Pointer reader = ReaderType::New();
        ReaderType::Pointer readerS = ReaderType::New();
        WriterType::Pointer writer = WriterType::New();
        //WriterType::Pointer writer_out = WriterType::New();
        stringstream ss;

        string name = "/home/gustavo/temp/cine_";
        string type = ".tif";

        ss<<name<<(i+1)<<type;

        string filename = ss.str();
        ss.str("");

        reader->SetFileName(filename);
        reader->Update();

        InternalImageType::Pointer val = reader->GetOutput();
        typedef itk::RescaleIntensityImageFilter<
                InternalImageType,
                OutputImageType >   CastFilterType;

        typedef   itk::CurvatureAnisotropicDiffusionImageFilter<
                InternalImageType,
                InternalImageType >  SmoothingFilterType;
        SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();

        typedef   itk::GradientMagnitudeRecursiveGaussianImageFilter<
                InternalImageType,
                InternalImageType >  GradientFilterType;
        typedef   itk::SigmoidImageFilter<
                InternalImageType,
                InternalImageType >  SigmoidFilterType;
        GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
        SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();

        sigmoid->SetOutputMinimum(sig_min);
        sigmoid->SetOutputMaximum(sig_max);

        typedef  itk::FastMarchingImageFilter<
                InternalImageType,
                InternalImageType >    FastMarchingFilterType;

        FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
        //const InternalImageType * inputImage = reader->GetOutput();
        const InternalImageType * inputImage = val;
        fastMarching->SetOutputRegion( inputImage->GetBufferedRegion() );
        fastMarching->SetOutputSpacing( inputImage->GetSpacing() );
        fastMarching->SetOutputOrigin( inputImage->GetOrigin() );
        fastMarching->SetOutputDirection( inputImage->GetDirection() );

        typedef  itk::GeodesicActiveContourLevelSetImageFilter< InternalImageType,
                InternalImageType >    GeodesicActiveContourFilterType;
        GeodesicActiveContourFilterType::Pointer geodesicActiveContour =
                GeodesicActiveContourFilterType::New();

        const double propagationScaling = propagation;
        //  Software Guide : BeginCodeSnippet
        geodesicActiveContour->SetPropagationScaling( propagationScaling );
        geodesicActiveContour->SetCurvatureScaling(curvature);
        geodesicActiveContour->SetAdvectionScaling(advection);

        geodesicActiveContour->SetMaximumRMSError(rms);
        geodesicActiveContour->SetNumberOfIterations(iterations);

        smoothing->SetInput(val);
        gradientMagnitude->SetInput( smoothing->GetOutput() );
        sigmoid->SetInput( gradientMagnitude->GetOutput() );
        geodesicActiveContour->SetInput(  fastMarching->GetOutput() );
        geodesicActiveContour->SetFeatureImage( sigmoid->GetOutput() );
        thresholder->SetInput( geodesicActiveContour->GetOutput() );
        //ImageOut = thresholder->GetOutput();
        //writer->SetInput( thresholder->GetOutput() );

        smoothing->SetTimeStep(timestep);
        smoothing->SetNumberOfIterations(it_dif);
        smoothing->SetConductanceParameter(conductance);

        const double sig = sigma;
        gradientMagnitude->SetSigma(sig);

        sigmoid->SetAlpha( alpha );
        sigmoid->SetBeta(  beta  );

        string line;
        string line2;
        int x;
        int y;
        int radius;
        if (myfile.is_open())
        {
            getline (myfile,line);
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
        if (myfile2.is_open())
        {
            getline(myfile2,line2);
            radius = atoi(line2.c_str());
        }

        stringstream segment;

        string nameS = "/home/gustavo/temp/segmented_";
        string typeS = ".tif";

        if(i<9)
            segment<<nameS<<"00"<<(i+1)<<typeS;
        if(i>=9 && i<99)
            segment<<nameS<<"0"<<(i+1)<<typeS;
        if(i>=99)
            segment<<nameS<<(i+1)<<typeS;
        string filenameS = segment.str();
        segment.str("");

        readerS->SetFileName(filenameS);
        readerS->Update();

        InternalImageType::Pointer valS = readerS->GetOutput();

        Utils utils;
        int seedX;
        int seedY;
        utils.GetSeed(valS, x, y, &seedX, &seedY);
        cout<<"SEED X: "<<seedX<<endl;
        cout<<"SEED Y: "<<seedY<<endl;

        typedef FastMarchingFilterType::NodeContainer  NodeContainer;
        typedef FastMarchingFilterType::NodeType       NodeType;
        NodeContainer::Pointer seeds = NodeContainer::New();
        InternalImageType::IndexType  seedPosition;
        seedPosition.SetElement(0, (seedX+4));
        seedPosition.SetElement(1,seedY);
        //seedPosition.SetElement(2,(int)z);
        cout<<"X, Y = "<<x<<" "<<y<<endl;
        NodeType node;
        const double seedValue = - distance;
        node.SetValue( seedValue );
        node.SetIndex( seedPosition );

        seeds->Initialize();
        seeds->InsertElement( 0, node );

        fastMarching->SetTrialPoints(  seeds  );

        fastMarching->SetSpeedConstant( 1.0 );

        CastFilterType::Pointer caster1 = CastFilterType::New();
        CastFilterType::Pointer caster2 = CastFilterType::New();
        CastFilterType::Pointer caster3 = CastFilterType::New();
        CastFilterType::Pointer caster4 = CastFilterType::New();
        WriterType::Pointer writer1 = WriterType::New();
        WriterType::Pointer writer2 = WriterType::New();
        WriterType::Pointer writer3 = WriterType::New();
        WriterType::Pointer writer4 = WriterType::New();
        caster1->SetInput( smoothing->GetOutput() );
        writer1->SetInput( caster1->GetOutput() );
        writer1->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput1_2.tif");
        caster1->SetOutputMinimum(   0 );
        caster1->SetOutputMaximum( 255 );
        writer1->Update();
        caster2->SetInput( gradientMagnitude->GetOutput() );
        writer2->SetInput( caster2->GetOutput() );
        writer2->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput2_2.tif");
        caster2->SetOutputMinimum(   0 );
        caster2->SetOutputMaximum( 255 );
        writer2->Update();
        caster3->SetInput( sigmoid->GetOutput() );
        writer3->SetInput( caster3->GetOutput() );
        writer3->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput3_2.tif");
        caster3->SetOutputMinimum(   0 );
        caster3->SetOutputMaximum( 255 );
        writer3->Update();
        caster4->SetInput( fastMarching->GetOutput() );
        writer4->SetInput( caster4->GetOutput() );
        writer4->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput4_2.tif");
        caster4->SetOutputMinimum(   0 );
        caster4->SetOutputMaximum( 255 );

        fastMarching->SetOutputSize(
                    reader->GetOutput()->GetBufferedRegion().GetSize() );
        reader->Update();

        std::cout << std::endl;
        std::cout << "Max. no. iterations: " << geodesicActiveContour->GetNumberOfIterations() << std::endl;
        std::cout << "Max. RMS error: " << geodesicActiveContour->GetMaximumRMSError() << std::endl;
        std::cout << std::endl;
        std::cout << "No. elpased iterations: " << geodesicActiveContour->GetElapsedIterations() << std::endl;
        std::cout << "RMS change: " << geodesicActiveContour->GetRMSChange() << std::endl;
        writer4->Update();

        typedef itk::ImageFileWriter< InternalImageType > InternalWriterType;
        InternalWriterType::Pointer mapWriter = InternalWriterType::New();
        mapWriter->SetInput( fastMarching->GetOutput() );
        mapWriter->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput4_2.mha");
        mapWriter->Update();
        InternalWriterType::Pointer speedWriter = InternalWriterType::New();
        speedWriter->SetInput( sigmoid->GetOutput() );
        speedWriter->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput3_2.mha");
        speedWriter->Update();
        InternalWriterType::Pointer gradientWriter = InternalWriterType::New();
        gradientWriter->SetInput( gradientMagnitude->GetOutput() );
        gradientWriter->SetFileName("/home/gustavo/temp/GeodesicActiveContourImageFilterOutput2_2.mha");
        gradientWriter->Update();

        CastFilterType::Pointer caster5 = CastFilterType::New();
        caster5->SetInput( reader->GetOutput() );
        caster5->SetOutputMinimum(   0 );
        caster5->SetOutputMaximum( 255 );
        caster5->Update();
        typedef itk::AndImageFilter <OutputImageType>  AndImageFilterType;
        AndImageFilterType::Pointer andFilter
                = AndImageFilterType::New();
        andFilter->SetInput(0, caster5->GetOutput());
        andFilter->SetInput(1, thresholder->GetOutput());
        andFilter->Update();

        stringstream ss2;

        string name2 = "/home/gustavo/temp/segmented2_";
        string type2 = ".tif";

        if(i<9)
            ss2<<name2<<"00"<<(i+1)<<type2;
        if(i>=9 && i<99)
            ss2<<name2<<"0"<<(i+1)<<type2;
        if(i>=99)
            ss2<<name2<<(i+1)<<type2;
        string filename2 = ss2.str();
        ss2.str("");

//        typedef itk::GradientMagnitudeImageFilter<OutputImageType, OutputImageType >  GradientType;
//        GradientType::Pointer gradient = GradientType::New();
//        gradient->SetInput(thresholder->GetOutput());
//        gradient->Update();

        try
        {
            writer->SetFileName(filename2);
            writer->SetInput(andFilter->GetOutput() );
            writer->Update();
        }
        catch( itk::ExceptionObject & excep )
        {
            std::cerr << "Exception caught !" << std::endl;
            std::cerr << excep << std::endl;
        }

    }
    myfile.close();

}
