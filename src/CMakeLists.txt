SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Werror")

SET (SPARK_SOURCES
     calc_rays.cpp
     characteristic.cpp
     grid.cpp
     main.cpp
     misc.cpp
     my_exceptions.cpp
     params.cpp
     read_params.cpp
    )

INCLUDE_DIRECTORIES (
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/misc
    $ENV{YAML_DIR}/include
    $ENV{BOOST_DIR}
    $ENV{GSL_DIR}/include
    $ENV{HDF_DIR}/include
    )

LINK_DIRECTORIES (
    $ENV{YAML_DIR}
    $ENV{HDF_DIR}/lib
    )

ADD_SUBDIRECTORY (utils)

ADD_EXECUTABLE (${PROJECT_NAME} ${SPARK_SOURCES})
TARGET_LINK_LIBRARIES (${PROJECT_NAME}
                        yaml-cpp
                        make_linelist
                        hdf5
                        hdf5_cpp
                      )
