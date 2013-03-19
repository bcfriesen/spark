AM_CPPFLAGS=${CPPFLAGS} -I../yaml-cpp/include
AM_LIBS = ../yaml-cpp/src/libyamlcpp.a
bin_PROGRAMS = spark
spark_SOURCES = calc_rays.cpp characteristic.cpp grid.cpp main.cpp misc.cpp my_exceptions.cpp params.cpp read_params.cpp
spark_LDADD = $(AM_LIBS)
