SET (RTCPP_SOURCES
    calc_rays.cpp
    charODE_dds.cpp
    charODE_dds_new.cpp
    charODE_ddr.cpp
    characteristic.cpp
    grid.cpp
    main.cpp
    misc.cpp
    my_exceptions.cpp
    )

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/src
                     $ENV{HOME}/yaml-cpp/include
                     $ENV{HOME}/odeint-v2
                     $ENV{HOME}/boost_1_51_0
                     )

LINK_DIRECTORIES ($ENV{HOME}/yaml-cpp/build)

ADD_EXECUTABLE (${PROJECT_NAME} ${RTCPP_SOURCES})
TARGET_LINK_LIBRARIES (${PROJECT_NAME} yaml-cpp)
