bin_PROGRAMS = spark
spark_CPPFLAGS=${CPPFLAGS} -I${top_srcdir}/yaml-cpp/include -I${BOOST_DIR}
spark_LIBS = ${top_builddir}/yaml-cpp/src/libyamlcpp.a
spark_SOURCES = \
				calc_rays.cpp \
				characteristic.cpp \
				grid.cpp \
				main.cpp \
				misc.cpp \
				my_exceptions.cpp \
				params.cpp \
				read_params.cpp
spark_LDADD = $(spark_LIBS)
