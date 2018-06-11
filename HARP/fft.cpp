#include "fft.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkWrapPadImageFilter.h"
#include "itkForwardFFTImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include <iostream>
#include "fft.h"
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

FFT::FFT()
{
}

using namespace std;

int Pad(int value){
    int coord = 0;
    if(value % 2 == 0){
        return coord;
    }
    else{
        coord = (value % 2);
        return coord;
    }
}

void FFT::Execute(ImageType::Pointer image){

    const unsigned int Dimension = 3;
    typedef float                                   FloatPixelType;
    typedef itk::Image< FloatPixelType, Dimension > FloatImageType;

    typedef unsigned char UnsignedCharPixelType;
    typedef itk::Image< UnsignedCharPixelType, Dimension > UnsignedCharImageType;
    // Some FFT filter implementations, like VNL's, need the image size to be a
    // multiple of small prime numbers.

    typedef itk::WrapPadImageFilter< FloatImageType, FloatImageType > PadFilterType;
    PadFilterType::Pointer padFilter = PadFilterType::New();
    padFilter->SetInput(image);

    PadFilterType::SizeType padding;
    const char *homedir;
    homedir = getpwuid(getuid())->pw_dir;
    string real;
    string imag;
    string mod;
    real = "/temp/fft/real.png";
    imag = "/temp/fft/imaginary.png";
    mod = "/temp/fft/modulus.png";
    string pathFileReal = homedir + real;
    string pathFileImag = homedir + imag;
    string pathFileMod = homedir + mod;
    // Input size is [48, 62, 42].  Pad to [48, 64, 48].


    int sizeX = image->GetLargestPossibleRegion().GetSize()[0];
    int sizeY = image->GetLargestPossibleRegion().GetSize()[1];
    int sizeZ = image->GetLargestPossibleRegion().GetSize()[2];

    cout<<"size x: "<<sizeX<<endl;
    cout<<"size y: "<<sizeY<<endl;
    cout<<"size z: "<<sizeZ<<endl;

    padding[0] = (-1)*Pad(sizeX);
    padding[1] = (-1)*Pad(sizeY);
    padding[2] = (-1)*Pad(sizeZ);
    padFilter->SetPadUpperBound( padding );
    padFilter->Update();

    typedef itk::ForwardFFTImageFilter< FloatImageType > FFTType;
    FFTType::Pointer fftFilter = FFTType::New();
    fftFilter->SetInput( padFilter->GetOutput() );
    typedef FFTType::OutputImageType FFTOutputImageType;
    // Extract the real part
    typedef itk::ComplexToRealImageFilter< FFTOutputImageType, FloatImageType> RealFilterType;
    RealFilterType::Pointer realFilter = RealFilterType::New();
    fftFilter->Update();
    realFilter->SetInput(fftFilter->GetOutput());
    typedef itk::RescaleIntensityImageFilter< FloatImageType, UnsignedCharImageType > RescaleFilterType;
    RescaleFilterType::Pointer realRescaleFilter = RescaleFilterType::New();
    realFilter->Update();
    realRescaleFilter->SetInput(realFilter->GetOutput());
    realRescaleFilter->SetOutputMinimum( itk::NumericTraits< UnsignedCharPixelType >::min() );
    realRescaleFilter->SetOutputMaximum( itk::NumericTraits< UnsignedCharPixelType >::max() );
    realRescaleFilter->Update();
    typedef itk::ImageFileWriter< UnsignedCharImageType > WriterType;
    WriterType::Pointer realWriter = WriterType::New();
    realWriter->SetFileName(pathFileReal.c_str());
    realWriter->SetInput( realRescaleFilter->GetOutput() );

    realWriter->Update();


    // Extract the imaginary part
    typedef itk::ComplexToImaginaryImageFilter< FFTOutputImageType, FloatImageType> ImaginaryFilterType;
    ImaginaryFilterType::Pointer imaginaryFilter = ImaginaryFilterType::New();
    imaginaryFilter->SetInput(fftFilter->GetOutput());
    RescaleFilterType::Pointer imaginaryRescaleFilter = RescaleFilterType::New();
    imaginaryFilter->Update();
    imaginaryRescaleFilter->SetInput(imaginaryFilter->GetOutput());
    imaginaryRescaleFilter->SetOutputMinimum( itk::NumericTraits< UnsignedCharPixelType >::min() );
    imaginaryRescaleFilter->SetOutputMaximum( itk::NumericTraits< UnsignedCharPixelType >::max() );
    imaginaryRescaleFilter->Update();
    WriterType::Pointer complexWriter = WriterType::New();
    complexWriter->SetFileName(pathFileImag.c_str());
    complexWriter->SetInput( imaginaryRescaleFilter->GetOutput() );

    complexWriter->Update();


    // Compute the magnitude
    typedef itk::ComplexToModulusImageFilter< FFTOutputImageType, FloatImageType> ModulusFilterType;
    ModulusFilterType::Pointer modulusFilter = ModulusFilterType::New();
    modulusFilter->SetInput(fftFilter->GetOutput());
    modulusFilter->Update();
    RescaleFilterType::Pointer magnitudeRescaleFilter = RescaleFilterType::New();
    magnitudeRescaleFilter->SetInput(modulusFilter->GetOutput());
    magnitudeRescaleFilter->SetOutputMinimum( itk::NumericTraits< UnsignedCharPixelType >::min() );
    magnitudeRescaleFilter->SetOutputMaximum( itk::NumericTraits< UnsignedCharPixelType >::max() );
    magnitudeRescaleFilter->Update();
    WriterType::Pointer magnitudeWriter = WriterType::New();
    magnitudeWriter->SetFileName(pathFileMod.c_str());
    magnitudeWriter->SetInput( magnitudeRescaleFilter->GetOutput() );

    magnitudeWriter->Update();

}
