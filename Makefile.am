# global compile and link options
AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(SIMD_FLAGS)
##AM_LDFLAGS = $(BOOST_REGEX_LDFLAGS) $(BOOST_REGEX_LIBS)
AM_LDFLAGS = -lboost_regex

# targets to build and install
lib_LTLIBRARIES = liblikely.la
bin_PROGRAMS = likelytest likelymc likelyinterp likelywsum likelyinteg likelyrand likelytricubic

# extra targets that should not be installed
noinst_PROGRAMS = demo1 demo2

# instructions for building the library
liblikely_la_SOURCES = \
	likely/FunctionMinimum.cc \
	likely/AbsEngine.cc \
	likely/EngineRegistry.cc \
	likely/MarkovChainEngine.cc \
	likely/Interpolator.cc \
	likely/Integrator.cc \
	likely/Random.cc \
	likely/AbsAccumulator.cc \
	likely/WeightedAccumulator.cc \
	likely/WeightedCombiner.cc \
	likely/QuantileAccumulator.cc \
	likely/TriCubicInterpolator.cc \
	likely/AbsBinning.cc \
	likely/test/TestLikelihood.cc

# library headers to install (nobase prefix preserves directories under bosslya)
# Anything that includes config.h should *not* be listed here.
# and is not part of the package public API.
nobase_include_HEADERS = \
	likely/likely.h \
	likely/types.h \
	likely/function.h \
	likely/function_impl.h \
	likely/RuntimeError.h \
	likely/FunctionMinimum.h \
	likely/AbsEngine.h \
	likely/EngineRegistry.h \
	likely/MarkovChainEngine.h \
	likely/Interpolator.h \
	likely/Integrator.h \
	likely/Random.h \
	likely/AbsAccumulator.h \
	likely/WeightedAccumulator.h \
	likely/WeightedCombiner.h \
	likely/QuantileAccumulator.h \
	likely/TriCubicInterpolator.h \
	likely/AbsBinning.h \
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

likelywsum_SOURCES = src/likelywsum.cc
likelywsum_DEPENDENCIES = $(lib_LIBRARIES)
likelywsum_LDADD = -lboost_program_options -L. -llikely

likelyinteg_SOURCES = src/likelyinteg.cc
likelyinteg_DEPENDENCIES = $(lib_LIBRARIES)
likelyinteg_LDADD = -L. -llikely

likelyrand_SOURCES = src/likelyrand.cc
likelyrand_DEPENDENCIES = $(lib_LIBRARIES)
likelyrand_LDADD = -L. -llikely

likelytricubic_SOURCES = src/likelytricubic.cc
likelytricubic_DEPENDENCIES = $(lib_LIBRARIES)
likelytricubic_LDADD = -lboost_program_options -L. -llikely

demo1_SOURCES = src/demo1.cc
demo1_DEPENDENCIES = $(lib_LIBRARIES)
demo1_LDADD = -L. -llikely

demo2_SOURCES = src/demo2.cc
demo2_DEPENDENCIES = $(lib_LIBRARIES)
demo2_LDADD = -L. -llikely
