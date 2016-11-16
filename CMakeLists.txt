cmake_minimum_required(VERSION 3.6)
project (gdal_GRIBAPI)

find_package(GDAL REQUIRED)
find_library(GRIBAPI_LIB grib_api)
find_file(GRIBAPI_INC grib_api.h)

include_directories(${GRIBAPI_INC})
include_directories(${GDAL_INLCUDE_DIR})
add_library(gdal_GRIBAPI MODULE gribdataset.cpp)
target_link_libraries(gdal_GRIBAPI ${GRIBAPI_LIB})

install (TARGETS gdal_GRIBAPI DESTINATION lib)
