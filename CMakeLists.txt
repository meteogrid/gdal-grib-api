cmake_minimum_required(VERSION 3.6)
project (gdal_GRIBAPI)

find_package(GDAL REQUIRED)
find_library(GRIBAPI_LIB grib_api)
find_file(GRIBAPI_HEADER grib_api.h)
get_filename_component(GRIBAPI_INCLUDE_DIR "$GRIBAPI_HEADER" DIRECTORY)

include_directories(${GRIBAPI_INCLUDE_DIR})
include_directories(${GDAL_INLCUDE_DIR})

add_library(gdal_GRIBAPI MODULE gribdataset.cpp)
set_target_properties(gdal_GRIBAPI PROPERTIES PREFIX "")
target_link_libraries(gdal_GRIBAPI ${GRIBAPI_LIB})

install (TARGETS gdal_GRIBAPI DESTINATION lib)
