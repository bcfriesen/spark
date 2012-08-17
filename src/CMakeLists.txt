SET (RTCPP_SOURCES
    main.cpp
    characteristic.cpp
    grid.cpp
    my_exceptions.cpp
    misc.cpp
    )

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src
    $ENV{HOME}/gsl/1.15/gnu/4.4.6/include)
LINK_DIRECTORIES($ENV{HOME}/gsl/1.15/gnu/4.4.6/lib)

SET (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

ADD_EXECUTABLE (${PROJECT_NAME} ${RTCPP_SOURCES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} gsl gslcblas)
