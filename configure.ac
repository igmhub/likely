#################################################################################
## Created 20-May-2011 by David Kirkby <dkirkby@uci.edu>
## Run autoreconf after making any change here.
#################################################################################

# This could probably be relaxed to an older version.
AC_PREREQ([2.68])

# Our package info.
AC_INIT([likely], [dev], [dkirkby@uci.edu])

# Checks for programs
AC_PROG_CXX

# Use SIMD compiler extensions when available
AX_EXT

# Initialize libtool, which adds --enable/disable-shared configure options.
# The boost.m4 macros used below also need this.
LT_INIT

# Checks for libraries (there is no portable way to check for C++ classes in
# the library, so we just check that the linker can find the library using 'main')
AC_CHECK_LIB([m],[cos],,
	AC_MSG_ERROR([Cannot find the math library.]))
AC_CHECK_LIB([lapack],[main],,
	AC_MSG_WARN([Cannot find the lapack library. Will try to continue without it.]))

# Following http://www.gentoo.org/proj/en/qa/automagic.xml below...

# Use 'configure --without-gsl' if GSL should not be used.
AC_ARG_WITH([gsl],
	AS_HELP_STRING([--without-gsl], [Build without the GSL library.]))
AS_IF([test "x$with_gsl" != "xno"], [
	AC_CHECK_LIB([gsl],[main],,
		AC_MSG_ERROR([Cannot find the gsl library.]))
])

# Use 'configure --without-minuit2' if Minuit2 should not be used.
AC_ARG_WITH([minuit2],
	AS_HELP_STRING([--without-minuit2], [Build without the Minuit2 library.]))
AS_IF([test "x$with_minuit2" != "xno"], [
	AC_CHECK_LIB([Minuit2], [main],,
		AC_MSG_ERROR([Cannot find the Minuit2 library.]))
])

# We need a recent version of boost
BOOST_REQUIRE([1.49])

# Required header-only boost packages
BOOST_BIND
BOOST_FOREACH
BOOST_FUNCTION
BOOST_UTILITY
BOOST_SMART_PTR
BOOST_TEST

# Required boost libraries
##BOOST_PROGRAM_OPTIONS
##BOOST_REGEX

# Configure automake
AC_CONFIG_FILES([Makefile])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
# See http://www.gnu.org/s/hello/manual/automake/maintainer_002dmode.html
AM_MAINTAINER_MODE
AC_CONFIG_HEADERS([config.h])

# Define automake variables that flag whether optional libraries should be used.
AM_CONDITIONAL([USE_GSL], [test "x$ac_cv_lib_gsl_main" = "xyes"])
AM_CONDITIONAL([USE_MINUIT2],[test "x$ac_cv_lib_Minuit2_main" = "xyes"])

AC_OUTPUT
