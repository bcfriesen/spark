SET (RTCPP_SOURCES
    main.cpp
    characteristic.cpp
    grid.cpp
    my_exceptions.cpp
    misc.cpp
    interpolate.cpp
    )

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)

SET (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

ADD_EXECUTABLE (${PROJECT_NAME} ${RTCPP_SOURCES})
