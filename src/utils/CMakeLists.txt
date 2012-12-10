SET(MAKE_LINELIST_SOURCES
    make_linelist_hdf.cpp
   )
INCLUDE_DIRECTORIES(
    ${PROJECT_SOURCE_DIR}/src/utils
    $ENV{BOOST_DIR}
    $ENV{HDF_DIR}/include
    )

ADD_LIBRARY(make_linelist
     ${MAKE_LINELIST_SOURCES}
    )
