#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([spark], [0.1], [friesen@ou.edu])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([src/config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

# Switch autoconf checks to C++ mode.
AC_LANG(C++)

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([cmath], [],
                 [AC_MSG_ERROR(cmath header required!)])
AC_CHECK_HEADERS([boost/numeric/odeint.hpp], [],
                 AC_MSG_ERROR(odeint.hpp required! (in Boost version >= 1.53.0)))

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

# DUMP
AC_CONFIG_FILES([
                 Makefile
                 src/Makefile
                 yaml-cpp/Makefile
                 yaml-cpp/src/Makefile])
AC_OUTPUT
