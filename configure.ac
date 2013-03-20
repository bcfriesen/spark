#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([spark], [0.1], [friesen@ou.edu])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([src/config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_RANLIB

# Switch autoconf checks to C++ mode.
AC_LANG([C++])

# Check for Boost
AC_ARG_WITH([boost],
            [AS_HELP_STRING([--with-boost=<dir>],
                            [Specify top-level directory containing Boost library.])])
AM_CONDITIONAL(USE_BOOST, [test ! "x$with_boost" = x])
AC_SUBST(BOOST_DIR, [$with_boost])

# Checks for libraries.

# Checks for cmath header
AC_CHECK_HEADERS([cmath], [],
                 [AC_MSG_ERROR(cmath header required!)])

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
