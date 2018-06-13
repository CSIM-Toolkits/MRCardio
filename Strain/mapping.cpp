#include "mapping.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkSymmetricForcesDemonsRegistrationFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"
#include <pwd.h>

using namespace std;

Mapping::Mapping()
{
    this->pw = getpwuid(getuid());

    this->homedir = this->pw->pw_dir;

    this->deformable = "/temp/strain/deformable_";
    this->vector = "/temp/strain/vector_";
    this->field = "/temp/strain/field_";
    this->pathDeformable = this->homedir + this->deformable;
    this->pathVector = this->homedir + this->vector;
    this->pathField = this->homedir + this->field;

}

/**
 * @brief Mapping::calcMapping
 * @param fixedImag
 * @param movingImag
 * @param index
 */
void Mapping::calcMapping(ImageType::Pointer fixedImag, ImageType::Pointer movingImag, int index){

    const unsigned int Dimension = 2;
    typedef unsigned short PixelType;
    string typeMha = ".mha";

    typedef itk::Image< PixelType, Dimension >  FixedImageType;
    typedef itk::Image< PixelType, Dimension >  MovingImageType;

    typedef float                                      InternalPixelType;
    typedef itk::Image< InternalPixelType, Dimension > InternalImageType;
    typedef itk::CastImageFilter< FixedImageType,
            InternalImageType >  FixedImageCasterType;
    typedef itk::CastImageFilter< MovingImageType,
            InternalImageType >  MovingImageCasterType;

    FixedImageCasterType::Pointer fixedImageCaster = FixedImageCasterType::New();
    MovingImageCasterType::Pointer movingImageCaster
            = MovingImageCasterType::New();

    fixedImageCaster->SetInput(fixedImag);
    movingImageCaster->SetInput(movingImag);

    typedef itk::HistogramMatchingImageFilter<
            InternalImageType,
            InternalImageType >   MatchingFilterType;
    MatchingFilterType::Pointer matcher = MatchingFilterType::New();

    matcher->SetInput( movingImageCaster->GetOutput() );
    matcher->SetReferenceImage( fixedImageCaster->GetOutput() );

    matcher->SetNumberOfHistogramLevels( 256 );
    matcher->SetNumberOfMatchPoints( 7 );

    matcher->ThresholdAtMeanIntensityOn();

    typedef itk::Vector< float, Dimension >           VectorPixelType;
    typedef itk::Image<  VectorPixelType, Dimension > DisplacementFieldType;
    typedef itk::SymmetricForcesDemonsRegistrationFilter<
            InternalImageType,
            InternalImageType,
            DisplacementFieldType> RegistrationFilterType;
    RegistrationFilterType::Pointer filter = RegistrationFilterType::New();

    filter->SetFixedImage( fixedImageCaster->GetOutput() );
    filter->SetMovingImage( matcher->GetOutput() );

    filter->SetNumberOfIterations( 5 );
    filter->SetStandardDeviations( 1.0 );

    filter->Update();

    typedef itk::WarpImageFilter<
            MovingImageType,
            MovingImageType,
            DisplacementFieldType  >     WarperType;
    typedef itk::LinearInterpolateImageFunction<
            MovingImageType,
            double          >  InterpolatorType;
    WarperType::Pointer warper = WarperType::New();
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    FixedImageType::Pointer fixedImage = fixedImag;

    warper->SetInput( movingImag );
    warper->SetInterpolator( interpolator );
    warper->SetOutputSpacing( fixedImage->GetSpacing() );
    warper->SetOutputOrigin( fixedImage->GetOrigin() );
    warper->SetOutputDirection( fixedImage->GetDirection() );

    warper->SetDisplacementField( filter->GetOutput() );

    typedef  unsigned char                           OutputPixelType;
    typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
    typedef itk::CastImageFilter<
            MovingImageType,
            OutputImageType >          CastFilterType;
    typedef itk::ImageFileWriter< OutputImageType >  WriterType;

    WriterType::Pointer      writer =  WriterType::New();
    CastFilterType::Pointer  caster =  CastFilterType::New();

    stringstream stringFileDeformable;
    stringFileDeformable<<pathDeformable<<(index+1)<<typeMha;

    string deformableFile = stringFileDeformable.str();
    stringFileDeformable.str("");

    writer->SetFileName(deformableFile);

    caster->SetInput( warper->GetOutput() );
    writer->SetInput( caster->GetOutput()   );
    writer->Update();

    typedef itk::ImageFileWriter< DisplacementFieldType > FieldWriterType;
    FieldWriterType::Pointer fieldWriter = FieldWriterType::New();

    stringstream stringFileField;
    stringFileField<<pathField<<(index+1)<<typeMha;

    string fieldFile = stringFileField.str();
    stringFileField.str("");

    fieldWriter->SetFileName(fieldFile);
    fieldWriter->SetInput( filter->GetOutput() );

    fieldWriter->Update();

    typedef DisplacementFieldType            VectorImage2DType;
    typedef DisplacementFieldType::PixelType Vector2DType;

    VectorImage2DType::ConstPointer vectorImage2D = filter->GetOutput();

    VectorImage2DType::RegionType  region2D = vectorImage2D->GetBufferedRegion();
    VectorImage2DType::IndexType   index2D  = region2D.GetIndex();
    VectorImage2DType::SizeType    size2D   = region2D.GetSize();


    typedef itk::Vector< float,       3 >  Vector3DType;
    typedef itk::Image< Vector3DType, 3 >  VectorImage3DType;

    typedef itk::ImageFileWriter< VectorImage3DType > VectorImage3DWriterType;

    VectorImage3DWriterType::Pointer writer3D = VectorImage3DWriterType::New();

    VectorImage3DType::Pointer vectorImage3D = VectorImage3DType::New();

    VectorImage3DType::RegionType  region3D;
    VectorImage3DType::IndexType   index3D;
    VectorImage3DType::SizeType    size3D;

    index3D[0] = index2D[0];
    index3D[1] = index2D[1];
    index3D[2] = 0;

    size3D[0]  = size2D[0];
    size3D[1]  = size2D[1];
    size3D[2]  = 1;

    region3D.SetSize( size3D );
    region3D.SetIndex( index3D );

    vectorImage3D->SetRegions( region3D );
    vectorImage3D->Allocate();

    typedef itk::ImageRegionConstIterator< VectorImage2DType > Iterator2DType;

    typedef itk::ImageRegionIterator< VectorImage3DType > Iterator3DType;

    Iterator2DType  it2( vectorImage2D, region2D );
    Iterator3DType  it3( vectorImage3D, region3D );

    it2.GoToBegin();
    it3.GoToBegin();

    Vector2DType vector2D;
    Vector3DType vector3D;

    vector3D[2] = 0; // set Z component to zero.

    while( !it2.IsAtEnd() ){
        vector2D = it2.Get();
        vector3D[0] = vector2D[0];
        vector3D[1] = vector2D[1];
        it3.Set( vector3D );
        ++it2;
        ++it3;
    }

    writer3D->SetInput( vectorImage3D );

    stringstream stringFileVector;
    stringFileVector<<pathVector<<(index+1)<<typeMha;

    string vectorFile = stringFileVector.str();
    stringFileVector.str("");

    writer3D->SetFileName(vectorFile);

    try
    {
        writer3D->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
        std::cerr << excp << std::endl;
    }

}
