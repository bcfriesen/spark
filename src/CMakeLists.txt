SET (SPARK_SOURCES
     calc_rays.cpp
     charODE_dds.cpp
     charODE_ddr.cpp
     characteristic.cpp
     grid.cpp
     main.cpp
     misc.cpp
     my_exceptions.cpp
    )

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/src
                     ${PROJECT_SOURCE_DIR}/src/EOS
                     $ENV{YAML_DIR}/include
                     $ENV{ODEINT_DIR}
                     $ENV{BOOST_DIR}
                     $ENV{GSL_DIR}/include
                     )

LINK_DIRECTORIES ($ENV{YAML_DIR})

ADD_SUBDIRECTORY (EOS)

ADD_EXECUTABLE (${PROJECT_NAME} ${SPARK_SOURCES})
TARGET_LINK_LIBRARIES (${PROJECT_NAME}
                        yaml-cpp
                        EOS
                      )
