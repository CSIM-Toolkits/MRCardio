#include "convert.h"
#include <string>
#include "itkImage.h"
#include <string.h>
#include <QString>
#include <iostream>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

using namespace std;
convert::convert()
{
    this->pw = getpwuid(getuid());

    this->homedir = this->pw->pw_dir;

    this->segmentedFinal = "/temp/cine_";

    this->pathSegmentedFinal = this->homedir + this->segmentedFinal;

    this->output = "/temp/output.tif";
    this->pathOutput = this->homedir + this->output;

}
typedef float PixelType;
typedef itk::Image<PixelType,3> ImageType;
ImageType::Pointer convert::ConvertImage(int first, int last){

    typedef float   PixelType;
      const unsigned int Dimension = 3;

      typedef itk::Image< PixelType, Dimension >  ImageType;
      typedef itk::ImageSeriesReader< ImageType > ReaderType;
      typedef itk::ImageFileWriter<   ImageType > WriterType;

      ReaderType::Pointer reader = ReaderType::New();
      WriterType::Pointer writer = WriterType::New();

      typedef itk::NumericSeriesFileNames    NameGeneratorType;

      NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
      string format = this->pathSegmentedFinal + "%03d.dcm";
      nameGenerator->SetSeriesFormat(format);

      nameGenerator->SetStartIndex( (first)+1 );
      nameGenerator->SetEndIndex( last );
      nameGenerator->SetIncrementIndex( 1 );
      std::vector<std::string> names = nameGenerator->GetFileNames();

      // List the files
      //
      std::vector<std::string>::iterator nit;
      for (nit = names.begin();
           nit != names.end();
           nit++)
        {
        std::cout << "File: " << (*nit).c_str() << std::endl;
        }

      reader->SetFileNames( names  );
      reader->Update();
      writer->SetFileName(this->pathOutput);
      writer->SetInput(reader->GetOutput());
      writer->Update();
      typedef itk::Image<float,3> ImageType2;
      ImageType2::Pointer imag_out = ImageType::New();
      imag_out = reader->GetOutput();
      return imag_out;


}
