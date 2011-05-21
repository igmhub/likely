# global compile and link options
AM_CPPFLAGS = -I@ROOTINCDIR@ $(BOOST_CPPFLAGS)
AM_LDFLAGS = -L@ROOTLIBDIR@

# targets to build and install
bin_PROGRAMS = likelytest
lib_LTLIBRARIES = liblikely.la

# instructions for building the library
liblikely_la_SOURCES = \
	likely/test/TestLikelihood.cc

# library headers to install (nobase prefix preserves directories under bosslya)
nobase_include_HEADERS = \
	likely/likely.h \
	likely/RuntimeError.h \
	likely/test/TestLikelihood.h

# instructions for building each program

likelytest_SOURCES = src/likelytest.cc
likelytest_DEPENDENCIES = $(lib_LIBRARIES)
likelytest_LDADD = @ROOTLIBS@ @ROOTAUXLIBS@ $(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
	$(BOOST_PROGRAM_OPTIONS_LIBS) -L. -llikely
