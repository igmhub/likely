# global compile and link options
AM_CPPFLAGS = $(BOOST_CPPFLAGS)
##AM_LDFLAGS = $(BOOST_REGEX_LDFLAGS) $(BOOST_REGEX_LIBS)
AM_LDFLAGS = -lboost_regex

# targets to build and install
lib_LTLIBRARIES = liblikely.la
bin_PROGRAMS = likelytest likelymc likelyinterp

# extra targets that should not be installed
noinst_PROGRAMS = demo1 demo2

# instructions for building the library
liblikely_la_SOURCES = \
	likely/FunctionMinimum.cc \
	likely/AbsEngine.cc \
	likely/EngineRegistry.cc \
	likely/MarkovChainEngine.cc \
	likely/Interpolator.cc \
	likely/Random.cc \
	likely/test/TestLikelihood.cc

# library headers to install (nobase prefix preserves directories under bosslya)
nobase_include_HEADERS = \
	likely/likely.h \
	likely/types.h \
	likely/RuntimeError.h \
	likely/FunctionMinimum.h \
	likely/AbsEngine.h \
	likely/EngineRegistry.h \
	likely/MarkovChainEngine.h \
	likely/Interpolator.h \
	likely/Random.h \
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
likelytest_LDADD = -lboost_program_options -L. -llikely

##	$(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) -L. -llikely

likelymc_SOURCES = src/likelymc.cc
likelymc_DEPENDENCIES = $(lib_LIBRARIES)
likelymc_LDADD = -lboost_program_options -L. -llikely

likelyinterp_SOURCES = src/likelyinterp.cc
likelyinterp_DEPENDENCIES = $(lib_LIBRARIES)
likelyinterp_LDADD = -lboost_program_options -L. -llikely

demo1_SOURCES = src/demo1.cc
demo1_DEPENDENCIES = $(lib_LIBRARIES)
demo1_LDADD = -L. -llikely

demo2_SOURCES = src/demo2.cc
demo2_DEPENDENCIES = $(lib_LIBRARIES)
demo2_LDADD = -L. -llikely
