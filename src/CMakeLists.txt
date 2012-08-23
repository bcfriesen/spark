SET (RTCPP_SOURCES
    main.cpp
    )

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/src
                     $ENV{HOME}/yaml-cpp/include
                     $ENV{HOME}/odeint-v2)

LINK_DIRECTORIES ($ENV{HOME}/yaml-cpp/build)

ADD_EXECUTABLE (${PROJECT_NAME} ${RTCPP_SOURCES})
TARGET_LINK_LIBRARIES (${PROJECT_NAME} yaml-cpp)
