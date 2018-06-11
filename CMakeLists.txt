cmake_minimum_required(VERSION 2.8.9)

project(CardioProg)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://www.example.com/Slicer/Extensions/CardioProg")
set(EXTENSION_CATEGORY "Registration")
set(EXTENSION_CONTRIBUTORS "Gustavo Canavaci Barizon (University of São Paulo)")
set(EXTENSION_DESCRIPTION "This is an extension for Segmentation and Registration of Cardiac Images")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/CardioProg.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/CardioProg/Screenshots/1.png")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(AD)
add_subdirectory(QEntropy)
add_subdirectory(Registration)
add_subdirectory(FractalDimension)

add_subdirectory(CardioSegmentation)
add_subdirectory(Filtering)
add_subdirectory(HoughTransform)
add_subdirectory(GeodesicActiveContour2D)
add_subdirectory(TextureFeatures)
add_subdirectory(HARP)
add_subdirectory(SampleEntropy)
add_subdirectory(Strain)
add_subdirectory(Tickness)
add_subdirectory(Extract)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
