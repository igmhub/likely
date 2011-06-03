# global compile and link options
AM_CPPFLAGS = $(BOOST_CPPFLAGS)

# targets to build and install
lib_LTLIBRARIES = liblikely.la

# extra targets that should not be installed
noinst_PROGRAMS = likelytest

# instructions for building the library
liblikely_la_SOURCES = \
	likely/FunctionMinimum.cc \
	likely/AbsEngine.cc \
	likely/MarkovChainEngine.cc \
	likely/test/TestLikelihood.cc

# library headers to install (nobase prefix preserves directories under bosslya)
nobase_include_HEADERS = \
	likely/likely.h \
	likely/types.h \
	likely/RuntimeError.h \
	likely/FunctionMinimum.h \
	likely/AbsEngine.h \
	likely/MarkovChainEngine.h \
	likely/test/TestLikelihood.h

# add GSL features when libgsl is available
if USE_GSL
liblikely_la_SOURCES += likely/GslEngine.cc likely/GslErrorHandler.cc
nobase_include_HEADERS += likely/GslEngine.h likely/GslErrorHandler.h
endif

# add Minuit2 features when libMinuit2 is available
if USE_MINUIT2
liblikely_la_SOURCES += likely/MinuitEngine.cc
nobase_include_HEADERS += likely/MinuitEngine.h
endif

# instructions for building each program

likelytest_SOURCES = src/likelytest.cc
likelytest_DEPENDENCIES = $(lib_LIBRARIES)
likelytest_LDADD = \
	$(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
	$(BOOST_REGEX_LDFLAGS) $(BOOST_REGEX_LIBS) \
	-L. -llikely
