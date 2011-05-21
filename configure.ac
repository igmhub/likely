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

# Initialize libtool, which adds --enable/disable-shared configure options.
# The boost.m4 macros used below also need this.
LT_INIT

# Checks for libraries.
##AC_CHECK_LIB([cfitsio],[fits_init_cfitsio],,
##	AC_MSG_ERROR([Cannot find the cfitsio library.]),[-lm])

# We need a recent version of ROOT
ROOT_PATH([5.28/00],,
	AC_MSG_ERROR([Cannot find a recent version of ROOT.]))

# We need a recent version of boost
BOOST_REQUIRE([1.45])

# Required header-only boost packages
BOOST_BIND
BOOST_FOREACH
BOOST_FUNCTION
BOOST_UTILITY
BOOST_SMART_PTR

# Required boost libraries
BOOST_PROGRAM_OPTIONS

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# Configure automake
AC_CONFIG_FILES([Makefile])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_OUTPUT
