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
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

#include "itkOtsuThresholdImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

using namespace std;

segment::segment()
{
    this->pw = getpwuid(getuid());

    this->homedir = this->pw->pw_dir;

    this->endocardium = "/temp/endocardium.txt";
    this->radius = "/temp/radius.txt";
    this->slices = "/temp/slices.txt";
    this->segmented = "/temp/segmented_";
    this->segmentedArea = "/temp/segmentedArea_";
    this->segmentedFinal = "/temp/segmentedFinal_";
    this->cine = "/temp/cine_";
    this->extractValues = "/temp/extractedValuesInternal.txt";
    this->extractValuesMyocardium = "/temp/extractedValuesMyocardium.txt";

    this->gacFilter = "/temp/GACFilter.tif";
    this->gacGradient = "/temp/GACGradient.tif";
    this->gacSigmoid = "/temp/GACSigmoid.tif";
    this->gacMap = "/temp/GACMap.tif";

    this->gacFilterFinal = "/temp/GACFilterFinal.tif";
    this->gacGradientFinal = "/temp/GACGradientFinal.tif";
    this->gacSigmoidFinal = "/temp/GACSigmoidFinal.tif";
    this->gacMapFinal = "/temp/GACMapFinal.tif";

    this->gacGradientMHA = "/temp/GACGradient.mha";
    this->gacSigmoidMHA = "/temp/GACSigmoid.mha";
    this->gacMapMHA = "/temp/GACMap.mha";

    this->gacGradientMHAFinal = "/temp/GACGradientFinal.mha";
    this->gacSigmoidMHAFinal = "/temp/GACSigmoidFinal.mha";
    this->gacMapMHAFinal = "/temp/GACMapFinal.mha";

    this->pathEndocardium = this->homedir + this->endocardium;
    this->pathRadius = this->homedir + this->radius;
    this->pathSlices = this->homedir + this->slices;
    this->pathSegmented = this->homedir + this->segmented;
    this->pathSegmentedArea = this->homedir + this->segmentedArea;
    this->pathSegmentedFinal = this->homedir + this->segmentedFinal;
    this->pathCine = this->homedir + this->cine;

    this->pathGacFilter = this->homedir + this->gacFilter;
    this->pathGacGradient = this->homedir + this->gacGradient;
    this->pathGacSigmoid = this->homedir + this->gacSigmoid;
    this->pathGacMap = this->homedir + this->gacMap;

    this->pathGacFilterFinal = this->homedir + this->gacFilterFinal;
    this->pathGacGradientFinal = this->homedir + this->gacGradientFinal;
    this->pathGacSigmoidFinal = this->homedir + this->gacSigmoidFinal;
    this->pathGacMapFinal = this->homedir + this->gacMapFinal;

    this->pathGacGradientMHA = this->homedir + this->gacGradientMHA;
    this->pathGacSigmoidMHA = this->homedir + this->gacSigmoidMHA;
    this->pathGacMapMHA = this->homedir + this->gacMapMHA;

    this->pathGacGradientMHAFinal = this->homedir + this->gacGradientMHAFinal;
    this->pathGacSigmoidMHAFinal = this->homedir + this->gacSigmoidMHAFinal;
    this->pathGacMapMHAFinal = this->homedir + this->gacMapMHAFinal;

    this->pathExtractValues = this->homedir + this->extractValues;
    this->pathExtractValuesMyocardium = this->homedir + this->extractValuesMyocardium;
}

void segment::InternalEC(int first,int last, double sigma, double sig_min, double sig_max, double propagation, double curvature, double advection, double rms, int iterations, double timestep, int it_dif, double conductance, double alpha, double beta, double distance){
    ifstream endocardiumFile (this->pathEndocardium.c_str());
    ofstream extractValues(this->pathExtractValues.c_str());
    if (extractValues.is_open()){
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
            WriterType::Pointer writerArea = WriterType::New();
            //WriterType::Pointer writer_out = WriterType::New();
            stringstream stringFileCine;

            string typeTiff = ".tif";

            stringFileCine<<this->pathCine<<(i+1)<<typeTiff;

            string cineFile = stringFileCine.str();
            stringFileCine.str("");

            reader->SetFileName(cineFile);
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
            gradientMagnitude->SetInput(smoothing->GetOutput());
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
            writer1->SetFileName(this->pathGacFilter);
            caster1->SetOutputMinimum(   0 );
            caster1->SetOutputMaximum( 255 );
            writer1->Update();
            caster2->SetInput( gradientMagnitude->GetOutput() );
            writer2->SetInput( caster2->GetOutput() );
            writer2->SetFileName(this->pathGacGradient);
            caster2->SetOutputMinimum(   0 );
            caster2->SetOutputMaximum( 255 );
            writer2->Update();
            caster3->SetInput( sigmoid->GetOutput() );
            writer3->SetInput( caster3->GetOutput() );
            writer3->SetFileName(this->pathGacSigmoid);
            caster3->SetOutputMinimum(   0 );
            caster3->SetOutputMaximum( 255 );
            writer3->Update();
            caster4->SetInput( fastMarching->GetOutput() );
            writer4->SetInput( caster4->GetOutput() );
            writer4->SetFileName(this->pathGacMap);
            caster4->SetOutputMinimum(   0 );
            caster4->SetOutputMaximum( 255 );

            fastMarching->SetOutputSize(
                        reader->GetOutput()->GetBufferedRegion().GetSize() );
            reader->Update();

            stringstream stringFileSegmented;

            if(i<9)
                stringFileSegmented<<this->pathSegmented<<"00"<<(i+1)<<typeTiff;
            if(i>=9 && i<99)
                stringFileSegmented<<this->pathSegmented<<"0"<<(i+1)<<typeTiff;
            if(i>=99)
                stringFileSegmented<<this->pathSegmented<<(i+1)<<typeTiff;
            string segmentedFile = stringFileSegmented.str();
            stringFileSegmented.str("");

            typedef itk::GradientMagnitudeImageFilter<OutputImageType, OutputImageType >  GradientType;
            GradientType::Pointer gradient = GradientType::New();
            gradient->SetInput(thresholder->GetOutput());
            gradient->Update();

            try
            {
                writer->SetFileName(segmentedFile);
                writer->SetInput( gradient->GetOutput() );
                writer->Update();
            }
            catch( itk::ExceptionObject & excep )
            {
                std::cerr << "Exception caught !" << std::endl;
                std::cerr << excep << std::endl;
            }

            stringstream stringFileAreaSegmented;

            if(i<9)
                stringFileAreaSegmented<<this->pathSegmentedArea<<"00"<<(i+1)<<typeTiff;
            if(i>=9 && i<99)
                stringFileAreaSegmented<<this->pathSegmentedArea<<"0"<<(i+1)<<typeTiff;
            if(i>=99)
                stringFileAreaSegmented<<this->pathSegmentedArea<<(i+1)<<typeTiff;
            string segmentedAreaFile = stringFileAreaSegmented.str();
            stringFileAreaSegmented.str("");

            try
            {
                writerArea->SetFileName(segmentedAreaFile);
                writerArea->SetInput( thresholder->GetOutput() );
                writerArea->Update();
            }
            catch( itk::ExceptionObject & excep )
            {
                std::cerr << "Exception caught !" << std::endl;
                std::cerr << excep << std::endl;
            }
            OutputImageType::Pointer imag_outArea = OutputImageType::New();
            imag_outArea = thresholder->GetOutput();
            OutputImageType::Pointer imag_outPerimeter = OutputImageType::New();
            imag_outPerimeter = gradient->GetOutput();
            Utils extract;
            double area = extract.GetArea(imag_outArea);
            double perimeter = extract.GetPerimeter(imag_outPerimeter);
            extractValues<<"Area of Image "<<i+1<<": "<<area<<endl;
            extractValues<<"Perimeter of Image "<<i+1<<": "<<perimeter<<endl;
            extractValues<<"----------------------"<<endl;

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
            mapWriter->SetFileName(this->pathGacMapMHA);
            mapWriter->Update();
            InternalWriterType::Pointer speedWriter = InternalWriterType::New();
            speedWriter->SetInput( sigmoid->GetOutput() );
            speedWriter->SetFileName(this->pathGacSigmoidMHA);
            speedWriter->Update();
            InternalWriterType::Pointer gradientWriter = InternalWriterType::New();
            gradientWriter->SetInput( gradientMagnitude->GetOutput() );
            gradientWriter->SetFileName(this->pathGacGradientMHA);
            gradientWriter->Update();
        }
    }
    extractValues.close();
    endocardiumFile.close();

}

void segment::MyocardiumEC(int first,int last, double sigma, double sig_min, double sig_max, double propagation, double curvature, double advection, double rms, int iterations, double timestep, int it_dif, double conductance, double alpha, double beta, double distance, bool up, bool down, bool left, bool hight){
    ifstream endocardiumFile (this->pathEndocardium.c_str());
    ofstream extractValuesMyocardium(this->pathExtractValuesMyocardium.c_str());
    ifstream radiusFile (this->pathRadius.c_str());
    if (extractValuesMyocardium.is_open()){
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
            WriterType::Pointer writerFillHole = WriterType::New();
            //WriterType::Pointer writer_out = WriterType::New();
            stringstream stringFileCine;

            string typeTiff = ".tif";

            stringFileCine<<this->pathCine<<(i+1)<<typeTiff;

            string cineFile = stringFileCine.str();
            stringFileCine.str("");

            reader->SetFileName(cineFile);
            reader->Update();

            InternalImageType::Pointer val = reader->GetOutput();

            typedef itk::OtsuThresholdImageFilter <InternalImageType, OutputImageType>
                    FilterType;
            FilterType::Pointer otsuFilter
                    = FilterType::New();
            otsuFilter->SetInput(val);
            otsuFilter->Update(); // To compute threshold

            typedef itk::VotingBinaryIterativeHoleFillingImageFilter<OutputImageType> FillHoleFilterType;
            FillHoleFilterType::InputSizeType radiusF;
            radiusF.Fill(3);
            FillHoleFilterType::Pointer filterFillHole = FillHoleFilterType::New();
            filterFillHole->SetInput( otsuFilter->GetOutput() );
            filterFillHole->SetRadius( radiusF );
            filterFillHole->SetMajorityThreshold(2);
            filterFillHole->SetBackgroundValue(255);
            filterFillHole->SetForegroundValue(0);
            filterFillHole->SetMaximumNumberOfIterations(10);

            stringstream stringFileSegmentedFillHole;

            if(i<9)
                stringFileSegmentedFillHole<<this->pathSegmentedFinal<<"FillHole"<<"00"<<(i+1)<<typeTiff;
            if(i>=9 && i<99)
                stringFileSegmentedFillHole<<this->pathSegmentedFinal<<"FillHole"<<"0"<<(i+1)<<typeTiff;
            if(i>=99)
                stringFileSegmentedFillHole<<this->pathSegmentedFinal<<"FillHole"<<(i+1)<<typeTiff;
            string fillHoleFile = stringFileSegmentedFillHole.str();
            stringFileSegmentedFillHole.str("");

            try
            {
                writerFillHole->SetFileName(fillHoleFile);
                writerFillHole->SetInput(filterFillHole->GetOutput() );
                writerFillHole->Update();
            }
            catch( itk::ExceptionObject & excep )
            {
                std::cerr << "Exception caught !" << std::endl;
                std::cerr << excep << std::endl;
            }

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
            gradientMagnitude->SetInput(smoothing->GetOutput());
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
            if (radiusFile.is_open())
            {
                getline(radiusFile,line2);
                radius = atoi(line2.c_str());
            }

            stringstream segment;

            if(i<9)
                segment<<this->pathSegmented<<"00"<<(i+1)<<typeTiff;
            if(i>=9 && i<99)
                segment<<this->pathSegmented<<"0"<<(i+1)<<typeTiff;
            if(i>=99)
                segment<<this->pathSegmented<<(i+1)<<typeTiff;
            string filenameSegmented = segment.str();
            segment.str("");

            readerS->SetFileName(filenameSegmented);
            readerS->Update();

            typedef itk::RescaleIntensityImageFilter<OutputImageType, InternalImageType >   FillCastType;

            FillCastType::Pointer filterFillCast = FillCastType::New();
            filterFillCast->SetInput( filterFillHole->GetOutput() );
            filterFillCast->SetOutputMinimum(   0 );
            filterFillCast->SetOutputMaximum( 255 );
            filterFillCast->Update();

            InternalImageType::Pointer valS = filterFillCast->GetOutput();

            int contSeed = 0;
            typedef FastMarchingFilterType::NodeContainer  NodeContainer;
            typedef FastMarchingFilterType::NodeType       NodeType;
            NodeContainer::Pointer seeds = NodeContainer::New();
            Utils utils;

            seeds->Initialize();

            if(up){
                int seedXUP;
                int seedYUP;
                utils.GetSeedUp(valS, x, y, &seedXUP, &seedYUP);
                InternalImageType::IndexType  seedPositionUp;
                seedPositionUp.SetElement(0,(seedXUP));
                seedPositionUp.SetElement(1,(seedYUP - 3));

                NodeType nodeUp;
                const double seedValue = - distance;
                nodeUp.SetValue( seedValue );
                nodeUp.SetIndex( seedPositionUp );

                seeds->InsertElement( contSeed, nodeUp );
                contSeed++;
            }
            if(down){
                int seedXDOWN;
                int seedYDOWN;
                utils.GetSeedDown(valS, x, y, &seedXDOWN, &seedYDOWN);
                InternalImageType::IndexType  seedPositionDown;
                seedPositionDown.SetElement(0,(seedXDOWN));
                seedPositionDown.SetElement(1,(seedYDOWN + 3));

                NodeType nodeDown;
                const double seedValue = - distance;
                nodeDown.SetValue( seedValue );
                nodeDown.SetIndex( seedPositionDown );

                seeds->InsertElement( contSeed, nodeDown );
                contSeed++;
            }
            if(left){
                int seedXLEFT;
                int seedYLEFT;
                utils.GetSeedLeft(valS, x, y, &seedXLEFT, &seedYLEFT);
                InternalImageType::IndexType  seedPositionLeft;
                seedPositionLeft.SetElement(0,(seedXLEFT - 3));
                seedPositionLeft.SetElement(1,(seedYLEFT));

                NodeType nodeLeft;
                const double seedValue = - distance;
                nodeLeft.SetValue( seedValue );
                nodeLeft.SetIndex( seedPositionLeft );

                seeds->InsertElement( contSeed, nodeLeft );
                contSeed++;
            }
            if(hight){
                int seedXHIGHT;
                int seedYHIGHT;
                utils.GetSeedHight(valS, x, y, &seedXHIGHT, &seedYHIGHT);
                InternalImageType::IndexType  seedPositionHight;
                seedPositionHight.SetElement(0,(seedXHIGHT + 3));
                seedPositionHight.SetElement(1,(seedYHIGHT));

                NodeType nodeHight;
                const double seedValue = - distance;
                nodeHight.SetValue( seedValue );
                nodeHight.SetIndex( seedPositionHight );

                seeds->InsertElement( contSeed, nodeHight );
                contSeed++;
            }

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
            caster1->SetInput(smoothing->GetOutput());
            writer1->SetInput( caster1->GetOutput() );
            writer1->SetFileName(this->pathGacFilterFinal);
            caster1->SetOutputMinimum(   0 );
            caster1->SetOutputMaximum( 255 );
            writer1->Update();
            caster2->SetInput( gradientMagnitude->GetOutput() );
            writer2->SetInput( caster2->GetOutput() );
            writer2->SetFileName(this->pathGacGradientFinal);
            caster2->SetOutputMinimum(   0 );
            caster2->SetOutputMaximum( 255 );
            writer2->Update();
            caster3->SetInput( sigmoid->GetOutput() );
            writer3->SetInput( caster3->GetOutput() );
            writer3->SetFileName(this->pathGacSigmoidFinal);
            caster3->SetOutputMinimum(   0 );
            caster3->SetOutputMaximum( 255 );
            writer3->Update();
            caster4->SetInput( fastMarching->GetOutput() );
            writer4->SetInput( caster4->GetOutput() );
            writer4->SetFileName(this->pathGacMapFinal);
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
            mapWriter->SetFileName(this->pathGacMapMHAFinal);
            mapWriter->Update();
            InternalWriterType::Pointer speedWriter = InternalWriterType::New();
            speedWriter->SetInput( sigmoid->GetOutput() );
            speedWriter->SetFileName(this->pathGacSigmoidMHAFinal);
            speedWriter->Update();
            InternalWriterType::Pointer gradientWriter = InternalWriterType::New();
            gradientWriter->SetInput( gradientMagnitude->GetOutput() );
            gradientWriter->SetFileName(this->pathGacGradientMHAFinal);
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

            stringstream stringFileSegmented;

            if(i<9)
                stringFileSegmented<<this->pathSegmentedFinal<<"00"<<(i+1)<<typeTiff;
            if(i>=9 && i<99)
                stringFileSegmented<<this->pathSegmentedFinal<<"0"<<(i+1)<<typeTiff;
            if(i>=99)
                stringFileSegmented<<this->pathSegmentedFinal<<(i+1)<<typeTiff;
            string segmentedFile = stringFileSegmented.str();
            stringFileSegmented.str("");

            typedef itk::GradientMagnitudeImageFilter<OutputImageType, OutputImageType >  GradientType;
            GradientType::Pointer gradient = GradientType::New();
            gradient->SetInput(thresholder->GetOutput());
            gradient->Update();

            try
            {
                writer->SetFileName(segmentedFile);
                writer->SetInput(andFilter->GetOutput() );
                writer->Update();
            }
            catch( itk::ExceptionObject & excep )
            {
                std::cerr << "Exception caught !" << std::endl;
                std::cerr << excep << std::endl;
            }

            stringstream stringFileAreaSegmented;

            if(i<9)
                stringFileAreaSegmented<<this->pathSegmentedArea<<"00"<<(i+1)<<typeTiff;
            if(i>=9 && i<99)
                stringFileAreaSegmented<<this->pathSegmentedArea<<"0"<<(i+1)<<typeTiff;
            if(i>=99)
                stringFileAreaSegmented<<this->pathSegmentedArea<<(i+1)<<typeTiff;
            string segmentedAreaFile = stringFileAreaSegmented.str();
            stringFileAreaSegmented.str("");

            WriterType::Pointer writerArea = WriterType::New();
            try
            {
                writerArea->SetFileName(segmentedAreaFile);
                writerArea->SetInput( thresholder->GetOutput() );
                writerArea->Update();
            }
            catch( itk::ExceptionObject & excep )
            {
                std::cerr << "Exception caught !" << std::endl;
                std::cerr << excep << std::endl;
            }

            OutputImageType::Pointer imag_outArea = OutputImageType::New();
            imag_outArea = thresholder->GetOutput();
            OutputImageType::Pointer imag_outPerimeter = OutputImageType::New();
            imag_outPerimeter = gradient->GetOutput();

            Utils extract;

            double area = extract.GetArea(imag_outArea);
            double perimeter = extract.GetPerimeter(imag_outPerimeter);
            extractValuesMyocardium<<"Area of Myocardium Image "<<i+1<<": "<<area<<endl;
            extractValuesMyocardium<<"Perimeter of Myocardium Image "<<i+1<<": "<<perimeter<<endl;
            extractValuesMyocardium<<"----------------------"<<endl;

        }
    }
    extractValuesMyocardium.close();
    endocardiumFile.close();

}

void segment::InternalELV(int first,int last, double sigma, double sig_min, double sig_max, double propagation, double curvature, double advection, double rms, int iterations, double timestep, int it_dif, double conductance, double alpha, double beta, double distance){
    ifstream endocardiumFile (this->pathEndocardium.c_str());
    ofstream extractValues(this->pathExtractValues.c_str());
    if (extractValues.is_open()){
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
            stringstream stringFileCine;

            string typeTiff = ".tif";

            stringFileCine<<this->pathCine<<(i+1)<<typeTiff;

            string cineFile = stringFileCine.str();
            stringFileCine.str("");

            reader->SetFileName(cineFile);
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
            gradientMagnitude->SetInput(smoothing->GetOutput());
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

            typedef FastMarchingFilterType::NodeContainer  NodeContainer;
            typedef FastMarchingFilterType::NodeType       NodeType;
            NodeContainer::Pointer seeds = NodeContainer::New();
            InternalImageType::IndexType  seedPosition;
            int seedX;
            int seedY;
            Utils utils;
            utils.GetCenter(val, &seedX, &seedY);
            seedPosition.SetElement(0,x);
            seedPosition.SetElement(1,y);
            //seedPosition.SetElement(2,(int)z);
            cout<<"X, Y = "<<seedX<<" "<<seedY<<endl;
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
            WriterType::Pointer writerArea = WriterType::New();
            caster1->SetInput(smoothing->GetOutput());
            writer1->SetInput( caster1->GetOutput() );
            writer1->SetFileName(this->pathGacFilter);
            caster1->SetOutputMinimum(   0 );
            caster1->SetOutputMaximum( 255 );
            writer1->Update();
            caster2->SetInput( gradientMagnitude->GetOutput() );
            writer2->SetInput( caster2->GetOutput() );
            writer2->SetFileName(this->pathGacGradient);
            caster2->SetOutputMinimum(   0 );
            caster2->SetOutputMaximum( 255 );
            writer2->Update();
            caster3->SetInput( sigmoid->GetOutput() );
            writer3->SetInput( caster3->GetOutput() );
            writer3->SetFileName(this->pathGacSigmoid);
            caster3->SetOutputMinimum(   0 );
            caster3->SetOutputMaximum( 255 );
            writer3->Update();
            caster4->SetInput( fastMarching->GetOutput() );
            writer4->SetInput( caster4->GetOutput() );
            writer4->SetFileName(this->pathGacMap);
            caster4->SetOutputMinimum(   0 );
            caster4->SetOutputMaximum( 255 );

            fastMarching->SetOutputSize(
                        reader->GetOutput()->GetBufferedRegion().GetSize() );
            reader->Update();

            stringstream stringFileSegmented;

            if(i<9)
                stringFileSegmented<<this->pathSegmented<<"00"<<(i+1)<<typeTiff;
            if(i>=9 && i<99)
                stringFileSegmented<<this->pathSegmented<<"0"<<(i+1)<<typeTiff;
            if(i>=99)
                stringFileSegmented<<this->pathSegmented<<(i+1)<<typeTiff;
            string segmentedFile = stringFileSegmented.str();
            stringFileSegmented.str("");

            typedef itk::GradientMagnitudeImageFilter<OutputImageType, OutputImageType >  GradientType;
            GradientType::Pointer gradient = GradientType::New();
            gradient->SetInput(thresholder->GetOutput());
            gradient->Update();

            try
            {
                writer->SetFileName(segmentedFile);
                writer->SetInput( gradient->GetOutput() );
                writer->Update();
            }
            catch( itk::ExceptionObject & excep )
            {
                std::cerr << "Exception caught !" << std::endl;
                std::cerr << excep << std::endl;
            }

            stringstream stringFileAreaSegmented;

            if(i<9)
                stringFileAreaSegmented<<this->pathSegmentedArea<<"00"<<(i+1)<<typeTiff;
            if(i>=9 && i<99)
                stringFileAreaSegmented<<this->pathSegmentedArea<<"0"<<(i+1)<<typeTiff;
            if(i>=99)
                stringFileAreaSegmented<<this->pathSegmentedArea<<(i+1)<<typeTiff;
            string segmentedAreaFile = stringFileAreaSegmented.str();
            stringFileAreaSegmented.str("");

            try
            {
                writerArea->SetFileName(segmentedAreaFile);
                writerArea->SetInput( thresholder->GetOutput() );
                writerArea->Update();
            }
            catch( itk::ExceptionObject & excep )
            {
                std::cerr << "Exception caught !" << std::endl;
                std::cerr << excep << std::endl;
            }
            OutputImageType::Pointer imag_outArea = OutputImageType::New();
            imag_outArea = thresholder->GetOutput();
            OutputImageType::Pointer imag_outPerimeter = OutputImageType::New();
            imag_outPerimeter = gradient->GetOutput();
            Utils extract;
            double area = extract.GetArea(imag_outArea);
            double perimeter = extract.GetPerimeter(imag_outPerimeter);
            extractValues<<"Area of Image "<<i+1<<": "<<area<<endl;
            extractValues<<"Perimeter of Image "<<i+1<<": "<<perimeter<<endl;
            extractValues<<"----------------------"<<endl;

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
            mapWriter->SetFileName(this->pathGacMapMHA);
            mapWriter->Update();
            InternalWriterType::Pointer speedWriter = InternalWriterType::New();
            speedWriter->SetInput( sigmoid->GetOutput() );
            speedWriter->SetFileName(this->pathGacSigmoidMHA);
            speedWriter->Update();
            InternalWriterType::Pointer gradientWriter = InternalWriterType::New();
            gradientWriter->SetInput( gradientMagnitude->GetOutput() );
            gradientWriter->SetFileName(this->pathGacGradientMHA);
            gradientWriter->Update();
        }
        extractValues.close();
        endocardiumFile.close();
    }
}

void segment::MyocardiumELV(int first,int last, double sigma, double sig_min, double sig_max, double propagation, double curvature, double advection, double rms, int iterations, double timestep, int it_dif, double conductance, double alpha, double beta, double distance, bool up, bool down, bool left, bool hight){
    ifstream endocardiumFile (this->pathEndocardium.c_str());
    ifstream radiusFile (this->pathRadius.c_str());
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
        stringstream stringFileCine;

        string typeTiff = ".tif";

        stringFileCine<<this->pathCine<<(i+1)<<typeTiff;

        string cineFile = stringFileCine.str();
        stringFileCine.str("");

        reader->SetFileName(cineFile);
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
        gradientMagnitude->SetInput(smoothing->GetOutput());
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
        if (radiusFile.is_open())
        {
            getline(radiusFile,line2);
            radius = atoi(line2.c_str());
        }

        stringstream segment;

        if(i<9)
            segment<<this->pathSegmented<<"00"<<(i+1)<<typeTiff;
        if(i>=9 && i<99)
            segment<<this->pathSegmented<<"0"<<(i+1)<<typeTiff;
        if(i>=99)
            segment<<this->pathSegmented<<(i+1)<<typeTiff;
        string filenameSegmented = segment.str();
        segment.str("");

        readerS->SetFileName(filenameSegmented);
        readerS->Update();

        InternalImageType::Pointer valS = readerS->GetOutput();

        int contSeed = 0;
        typedef FastMarchingFilterType::NodeContainer  NodeContainer;
        typedef FastMarchingFilterType::NodeType       NodeType;
        NodeContainer::Pointer seeds = NodeContainer::New();
        Utils utils;

        seeds->Initialize();

        if(up){
            int seedXUP;
            int seedYUP;
            utils.GetSeedUp(valS, x, y, &seedXUP, &seedYUP);
            InternalImageType::IndexType  seedPositionUp;
            seedPositionUp.SetElement(0,(seedXUP));
            seedPositionUp.SetElement(1,(seedYUP - 3));

            NodeType nodeUp;
            const double seedValue = - distance;
            nodeUp.SetValue( seedValue );
            nodeUp.SetIndex( seedPositionUp );

            seeds->InsertElement( contSeed, nodeUp );
            contSeed++;
        }
        if(down){
            int seedXDOWN;
            int seedYDOWN;
            utils.GetSeedDown(valS, x, y, &seedXDOWN, &seedYDOWN);
            InternalImageType::IndexType  seedPositionDown;
            seedPositionDown.SetElement(0,(seedXDOWN));
            seedPositionDown.SetElement(1,(seedYDOWN + 3));

            NodeType nodeDown;
            const double seedValue = - distance;
            nodeDown.SetValue( seedValue );
            nodeDown.SetIndex( seedPositionDown );

            seeds->InsertElement( contSeed, nodeDown );
            contSeed++;
        }
        if(left){
            int seedXLEFT;
            int seedYLEFT;
            utils.GetSeedLeft(valS, x, y, &seedXLEFT, &seedYLEFT);
            InternalImageType::IndexType  seedPositionLeft;
            seedPositionLeft.SetElement(0,(seedXLEFT - 3));
            seedPositionLeft.SetElement(1,(seedYLEFT));

            NodeType nodeLeft;
            const double seedValue = - distance;
            nodeLeft.SetValue( seedValue );
            nodeLeft.SetIndex( seedPositionLeft );

            seeds->InsertElement( contSeed, nodeLeft );
            contSeed++;
        }
        if(hight){
            int seedXHIGHT;
            int seedYHIGHT;
            utils.GetSeedHight(valS, x, y, &seedXHIGHT, &seedYHIGHT);
            InternalImageType::IndexType  seedPositionHight;
            seedPositionHight.SetElement(0,(seedXHIGHT + 3));
            seedPositionHight.SetElement(1,(seedYHIGHT));

            NodeType nodeHight;
            const double seedValue = - distance;
            nodeHight.SetValue( seedValue );
            nodeHight.SetIndex( seedPositionHight );

            seeds->InsertElement( contSeed, nodeHight );
            contSeed++;
        }
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
        caster1->SetInput(smoothing->GetOutput());
        writer1->SetInput( caster1->GetOutput() );
        writer1->SetFileName(this->pathGacFilterFinal);
        caster1->SetOutputMinimum(   0 );
        caster1->SetOutputMaximum( 255 );
        writer1->Update();
        caster2->SetInput( gradientMagnitude->GetOutput() );
        writer2->SetInput( caster2->GetOutput() );
        writer2->SetFileName(this->pathGacGradientFinal);
        caster2->SetOutputMinimum(   0 );
        caster2->SetOutputMaximum( 255 );
        writer2->Update();
        caster3->SetInput( sigmoid->GetOutput() );
        writer3->SetInput( caster3->GetOutput() );
        writer3->SetFileName(this->pathGacSigmoidFinal);
        caster3->SetOutputMinimum(   0 );
        caster3->SetOutputMaximum( 255 );
        writer3->Update();
        caster4->SetInput( fastMarching->GetOutput() );
        writer4->SetInput( caster4->GetOutput() );
        writer4->SetFileName(this->pathGacMapFinal);
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
        mapWriter->SetFileName(this->pathGacMapMHAFinal);
        mapWriter->Update();
        InternalWriterType::Pointer speedWriter = InternalWriterType::New();
        speedWriter->SetInput( sigmoid->GetOutput() );
        speedWriter->SetFileName(this->pathGacSigmoidMHAFinal);
        speedWriter->Update();
        InternalWriterType::Pointer gradientWriter = InternalWriterType::New();
        gradientWriter->SetInput( gradientMagnitude->GetOutput() );
        gradientWriter->SetFileName(this->pathGacGradientMHAFinal);
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

        stringstream stringFileSegmented;

        if(i<9)
            stringFileSegmented<<this->pathSegmentedFinal<<"00"<<(i+1)<<typeTiff;
        if(i>=9 && i<99)
            stringFileSegmented<<this->pathSegmentedFinal<<"0"<<(i+1)<<typeTiff;
        if(i>=99)
            stringFileSegmented<<this->pathSegmentedFinal<<(i+1)<<typeTiff;
        string segmentedFile = stringFileSegmented.str();
        stringFileSegmented.str("");

        //        typedef itk::GradientMagnitudeImageFilter<OutputImageType, OutputImageType >  GradientType;
        //        GradientType::Pointer gradient = GradientType::New();
        //        gradient->SetInput(thresholder->GetOutput());
        //        gradient->Update();

        try
        {
            writer->SetFileName(segmentedFile);
            writer->SetInput(andFilter->GetOutput() );
            writer->Update();
        }
        catch( itk::ExceptionObject & excep )
        {
            std::cerr << "Exception caught !" << std::endl;
            std::cerr << excep << std::endl;
        }

    }
    endocardiumFile.close();

}
