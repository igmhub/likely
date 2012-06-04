// Created 08-Jun-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/Random.h"
#include "likely/RuntimeError.h"

#include "boost/random/uniform_01.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/random/uniform_int_distribution.hpp"
#include "boost/random/variate_generator.hpp"
#include "boost/lexical_cast.hpp"

#include <cmath>
#include <cstdlib>
#include <string>

#include "config.h" // defines HAVE_SSE2 when appropriate, thanks to the AX_EXT m4 macro
#define MEXP 19937
//#define MEXP 11213 // 7% faster in the non-SIMD version, but shorter period
#include "SFMT/SFMT.c"

namespace local = likely;

local::Random::Random() :
_uniform(boost::variate_generator<boost::mt19937&, boost::uniform_01<> >
    (_generator, boost::uniform_01<>())),
_gauss(boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
    (_generator, boost::normal_distribution<>(0,1)))
{
}

local::Random &local::Random::instance() {
    static Random *_instance = new Random();
    return *_instance;
}

void local::Random::setSeed(int seedValue) {
    _generator.seed(seedValue);
    init_gen_rand(seedValue);
}

int local::Random::getInteger(int min, int max) {
    boost::random::uniform_int_distribution<> dist(min,max);
    return dist(_generator);
}

void local::Random::partialShuffle(std::vector<int> &sample, int size) {
    int ssize(sample.size());
    if(size <= 0 || size > ssize) {
        throw RuntimeError("Random::partialShuffle: expected 0 < size <= sample.size().");
    }
    for(int pos = 0; pos < size; ++pos) {
        int pick = getInteger(pos+1,ssize-1);
        std::swap(sample[pos],sample[pick]);
    }
}

void local::Random::sampleWithReplacement(std::vector<int> &sample, int size) {
    if(size <= 0) {
        throw RuntimeError("Random::sampleWithReplacement: expected size > 0.");
    }
    int ssize(sample.size());
    sample.assign(ssize,0);
    boost::random::uniform_int_distribution<> dist(0,ssize-1);
    for(int trial = 0; trial < size; ++trial) {
        sample[dist(_generator)]++;
    }
}

float local::Random::getFastUniform() {
    return genrand_res53();
}

void *local::allocateAlignedArray(std::size_t byteSize) {
    void *array(0);
#if defined(__APPLE__) || \
    (defined(__FreeBSD__) && __FreeBSD__ >= 3 && __FreeBSD__ <= 6)
    array = malloc(byteSize);
    if (array == NULL) {
        throw RuntimeError("allocateAlignedArray: Apple malloc failed.");
    }
#elif defined(_POSIX_C_SOURCE)
    if (posix_memalign((void **)&array, 16, byteSize) != 0) {
        throw RuntimeError("allocateAlignedArray: posix_memalign failed.");
    }
#elif defined(__GNUC__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 3))
    array = (uint32_t*)memalign(16, byteSize);
    if (array == NULL) {
        throw RuntimeError("allocateAlignedArray: GNUC memalign failed.");
    }
#else /* in this case, gcc doesn't support SSE2 */
    array = malloc(byteSize);
    if (array == NULL) {
        throw RuntimeError("allocateAlignedArray: default malloc failed.");
    }
#endif
    return array;
}

boost::shared_array<float> local::allocateAlignedFloatArray(std::size_t size) {
    assert(sizeof(float) == sizeof(uint32_t));
    float *fbuffer = (float*)allocateAlignedArray(size*sizeof(uint32_t));
    return boost::shared_array<float>(fbuffer,std::ptr_fun(free));
}

boost::shared_array<double> local::allocateAlignedDoubleArray(std::size_t size) {
    assert(sizeof(double) == sizeof(uint64_t));
    double *dbuffer = (double*)allocateAlignedArray(size*sizeof(uint64_t));
    return boost::shared_array<double>(dbuffer,std::ptr_fun(free));
}

std::size_t local::Random::_initializeFill(std::size_t nrandom,
int seed, int stride, int minimum) {
    if(nrandom <= 0) {
        throw RuntimeError("Random: expected nrandom > 0.");
    }
    // Increase the number generated?
    std::size_t ngen(nrandom);
    if(ngen < minimum) {
        // Round up to the minimum size.
        ngen = minimum;
    }
    else if(ngen % stride) {
        // Round up for alignment.
        ngen += stride - (ngen % stride);
    }
    // Set the random seed.
    init_gen_rand(seed);
    if(!initialized || idx != N32) {
        throw RuntimeError("Random: init_gen_rand failed.");
    }
    return ngen;
}

boost::shared_array<double> local::Random::fillDoubleArrayUniform(std::size_t &nrandom, int seed) {
    // Get the number of random 64-bit integers to generate.
    nrandom = _initializeFill(nrandom,seed,2,N64);
    // Allocate the shared array.
    boost::shared_array<double> sarray = allocateAlignedDoubleArray(nrandom);
    double *array = sarray.get();
    // Fill the array with random bits.
    gen_rand_array((w128_t *)array, nrandom/2);
    idx = N32;
#if defined(BIG_ENDIAN64)
    swap((w128_t *)array, nrandom /2);
#endif
    uint64_t *ptr((uint64_t*)array);
    for(int i = 0; i < nrandom; ++i) {
        array[i] = to_res53(*ptr++);
    }
    return sarray;
}

boost::shared_array<double> local::Random::fillDoubleArrayNormal(std::size_t &nrandom, int seed) {
    // Round nrandom up to an even number to simplify alignment issues.
    if(nrandom % 2) nrandom++;
    // Get the number of random 32-bit integers to generate.
    std::size_t ngen = _initializeFill(nrandom,seed,4,N32);
    // Will this fit within an array of nrandom doubles? Ensure that nrandom >= ngen/2.
    if(ngen > 2*nrandom) nrandom = ngen/2;
    // Allocate the shared array with enough space for ngen quad-aligned 32-bit random integers
    // and nrandom 64-bit doubles.
    boost::shared_array<double> sarray = allocateAlignedDoubleArray(nrandom);
    double *array = sarray.get();
    // Calculate the 64-bit offset for filling the array in the top of the output array.
    int offset = nrandom - ngen/2;
    // Fill the array with random bits.
    gen_rand_array((w128_t *)(array+offset), ngen/4);
    idx = N32;    
    // Calculate where to start reading the 32-bit random integers so we will not
    // overwrite them as we save the new double values. Step n involves reading the next
    // 32-bit int from [offset+n] and writing the new double into [2n] and [2n+1], so
    // we require that the next 32-bit int not be clobbered, i.e., offset+n+1 > 2n+1,
    // or simply offset > n for all n. Since n covers the range 0..(nrandom-1), take
    // offset = nrandom. This will always fit since it is a 32-bit offset and the array
    // has space for at least nrandom 64-bit doubles.
    uint32_t *ptr((uint32_t*)array+nrandom);
    // Read random integers and convert them to normally distributed doubles.
    for(int index = 0; index < nrandom; ++index) {
        array[index] = _zigguratConvert(*ptr++);
    }
    return sarray;
}

boost::shared_array<float> local::Random::fillFloatArrayNormal(std::size_t &nrandom, int seed) {
    // Get the number of random 32-bit integers to generate.
    nrandom = _initializeFill(nrandom,seed,4,N32);
    // Allocate the shared array
    boost::shared_array<float> sarray = allocateAlignedFloatArray(nrandom);
    float *array = sarray.get();
    // Fill the array with random bits.
    gen_rand_array((w128_t *)array, nrandom/4);
    idx = N32;
    // Read random integers and convert them to normally distributed floats.
    uint32_t *ptr((uint32_t*)array);
    for(int index = 0; index < nrandom; ++index) {
        array[index] = (float)_zigguratConvert(*ptr++);
    }
    return sarray;
}

/* position of right-most step */
#define PARAM_R 3.44428647676

double local::Random::_zigguratConvert(uint32_t U) {
    uint32_t i,j,sign;
    double x, y;
    while(1) {
        i = U & 0x0000007F;		/* 7 bit to choose the step */
        sign = U & 0x00000080;	/* 1 bit for the sign */
        j = U>>8;   		    /* 24 bit for the x-value */

        x = j*_ziggurat_wtab[i];
        if (j < _ziggurat_ktab[i])  break;

        if (i<127) {
            double  y0, y1;
            y0 = _ziggurat_ytab[i];
            y1 = _ziggurat_ytab[i+1];
            y = y1+(y0-y1)*genrand_res53();
        }
        else {
            x = PARAM_R - std::log(1.0-genrand_res53())/PARAM_R;
            y = std::exp(-PARAM_R*(x-0.5*PARAM_R))*genrand_res53();
        }
        if (y < std::exp(-0.5*x*x))  break;
        // If we get here, we need a new 32-bit random number in U.
        // We actually generate a 64-bit random integer to stay in synch.
        U = gen_rand64() & 0xffffffff;
    }
    return sign ? +x : -x;
}

/* tabulated values for the heigt of the Ziggurat levels */
const double local::Random::_ziggurat_ytab[128] = {
  1, 0.963598623011, 0.936280813353, 0.913041104253,
  0.892278506696, 0.873239356919, 0.855496407634, 0.838778928349,
  0.822902083699, 0.807732738234, 0.793171045519, 0.779139726505,
  0.765577436082, 0.752434456248, 0.739669787677, 0.727249120285,
  0.715143377413, 0.703327646455, 0.691780377035, 0.68048276891,
  0.669418297233, 0.65857233912, 0.647931876189, 0.637485254896,
  0.62722199145, 0.617132611532, 0.607208517467, 0.597441877296,
  0.587825531465, 0.578352913803, 0.569017984198, 0.559815170911,
  0.550739320877, 0.541785656682, 0.532949739145, 0.524227434628,
  0.515614886373, 0.507108489253, 0.498704867478, 0.490400854812,
  0.482193476986, 0.47407993601, 0.466057596125, 0.458123971214,
  0.450276713467, 0.442513603171, 0.434832539473, 0.427231532022,
  0.419708693379, 0.41226223212, 0.404890446548, 0.397591718955,
  0.390364510382, 0.383207355816, 0.376118859788, 0.369097692334,
  0.362142585282, 0.355252328834, 0.348425768415, 0.341661801776,
  0.334959376311, 0.328317486588, 0.321735172063, 0.31521151497,
  0.308745638367, 0.302336704338, 0.29598391232, 0.289686497571,
  0.283443729739, 0.27725491156, 0.271119377649, 0.265036493387,
  0.259005653912, 0.253026283183, 0.247097833139, 0.241219782932,
  0.235391638239, 0.229612930649, 0.223883217122, 0.218202079518,
  0.212569124201, 0.206983981709, 0.201446306496, 0.195955776745,
  0.190512094256, 0.185114984406, 0.179764196185, 0.174459502324,
  0.169200699492, 0.1639876086, 0.158820075195, 0.153697969964,
  0.148621189348, 0.143589656295, 0.138603321143, 0.133662162669,
  0.128766189309, 0.123915440582, 0.119109988745, 0.114349940703,
  0.10963544023, 0.104966670533, 0.100343857232, 0.0957672718266,
  0.0912372357329, 0.0867541250127, 0.082318375932, 0.0779304915295,
  0.0735910494266, 0.0693007111742, 0.065060233529, 0.0608704821745,
  0.056732448584, 0.05264727098, 0.0486162607163, 0.0446409359769,
  0.0407230655415, 0.0368647267386, 0.0330683839378, 0.0293369977411,
  0.0256741818288, 0.0220844372634, 0.0185735200577, 0.0151490552854,
  0.0118216532614, 0.00860719483079, 0.00553245272614, 0.00265435214565
};

/* tabulated values for 2^24 times x[i]/x[i+1],
 * used to accept for U*x[i+1]<=x[i] without any floating point operations */
const uint32_t local::Random::_ziggurat_ktab[128] = {
  0, 12590644, 14272653, 14988939,
  15384584, 15635009, 15807561, 15933577,
  16029594, 16105155, 16166147, 16216399,
  16258508, 16294295, 16325078, 16351831,
  16375291, 16396026, 16414479, 16431002,
  16445880, 16459343, 16471578, 16482744,
  16492970, 16502368, 16511031, 16519039,
  16526459, 16533352, 16539769, 16545755,
  16551348, 16556584, 16561493, 16566101,
  16570433, 16574511, 16578353, 16581977,
  16585398, 16588629, 16591685, 16594575,
  16597311, 16599901, 16602354, 16604679,
  16606881, 16608968, 16610945, 16612818,
  16614592, 16616272, 16617861, 16619363,
  16620782, 16622121, 16623383, 16624570,
  16625685, 16626730, 16627708, 16628619,
  16629465, 16630248, 16630969, 16631628,
  16632228, 16632768, 16633248, 16633671,
  16634034, 16634340, 16634586, 16634774,
  16634903, 16634972, 16634980, 16634926,
  16634810, 16634628, 16634381, 16634066,
  16633680, 16633222, 16632688, 16632075,
  16631380, 16630598, 16629726, 16628757,
  16627686, 16626507, 16625212, 16623794,
  16622243, 16620548, 16618698, 16616679,
  16614476, 16612071, 16609444, 16606571,
  16603425, 16599973, 16596178, 16591995,
  16587369, 16582237, 16576520, 16570120,
  16562917, 16554758, 16545450, 16534739,
  16522287, 16507638, 16490152, 16468907,
  16442518, 16408804, 16364095, 16301683,
  16207738, 16047994, 15704248, 15472926
};

/* tabulated values of 2^{-24}*x[i] */
const double local::Random::_ziggurat_wtab[128] = {
  1.62318314817e-08, 2.16291505214e-08, 2.54246305087e-08, 2.84579525938e-08,
  3.10340022482e-08, 3.33011726243e-08, 3.53439060345e-08, 3.72152672658e-08,
  3.8950989572e-08, 4.05763964764e-08, 4.21101548915e-08, 4.35664624904e-08,
  4.49563968336e-08, 4.62887864029e-08, 4.75707945735e-08, 4.88083237257e-08,
  5.00063025384e-08, 5.11688950428e-08, 5.22996558616e-08, 5.34016475624e-08,
  5.44775307871e-08, 5.55296344581e-08, 5.65600111659e-08, 5.75704813695e-08,
  5.85626690412e-08, 5.95380306862e-08, 6.04978791776e-08, 6.14434034901e-08,
  6.23756851626e-08, 6.32957121259e-08, 6.42043903937e-08, 6.51025540077e-08,
  6.59909735447e-08, 6.68703634341e-08, 6.77413882848e-08, 6.8604668381e-08,
  6.94607844804e-08, 7.03102820203e-08, 7.11536748229e-08, 7.1991448372e-08,
  7.2824062723e-08, 7.36519550992e-08, 7.44755422158e-08, 7.52952223703e-08,
  7.61113773308e-08, 7.69243740467e-08, 7.77345662086e-08, 7.85422956743e-08,
  7.93478937793e-08, 8.01516825471e-08, 8.09539758128e-08, 8.17550802699e-08,
  8.25552964535e-08, 8.33549196661e-08, 8.41542408569e-08, 8.49535474601e-08,
  8.57531242006e-08, 8.65532538723e-08, 8.73542180955e-08, 8.8156298059e-08,
  8.89597752521e-08, 8.97649321908e-08, 9.05720531451e-08, 9.138142487e-08,
  9.21933373471e-08, 9.30080845407e-08, 9.38259651738e-08, 9.46472835298e-08,
  9.54723502847e-08, 9.63014833769e-08, 9.71350089201e-08, 9.79732621669e-08,
  9.88165885297e-08, 9.96653446693e-08, 1.00519899658e-07, 1.0138063623e-07,
  1.02247952126e-07, 1.03122261554e-07, 1.04003996769e-07, 1.04893609795e-07,
  1.05791574313e-07, 1.06698387725e-07, 1.07614573423e-07, 1.08540683296e-07,
  1.09477300508e-07, 1.1042504257e-07, 1.11384564771e-07, 1.12356564007e-07,
  1.13341783071e-07, 1.14341015475e-07, 1.15355110887e-07, 1.16384981291e-07,
  1.17431607977e-07, 1.18496049514e-07, 1.19579450872e-07, 1.20683053909e-07,
  1.21808209468e-07, 1.2295639141e-07, 1.24129212952e-07, 1.25328445797e-07,
  1.26556042658e-07, 1.27814163916e-07, 1.29105209375e-07, 1.30431856341e-07,
  1.31797105598e-07, 1.3320433736e-07, 1.34657379914e-07, 1.36160594606e-07,
  1.37718982103e-07, 1.39338316679e-07, 1.41025317971e-07, 1.42787873535e-07,
  1.44635331499e-07, 1.4657889173e-07, 1.48632138436e-07, 1.50811780719e-07,
  1.53138707402e-07, 1.55639532047e-07, 1.58348931426e-07, 1.61313325908e-07,
  1.64596952856e-07, 1.68292495203e-07, 1.72541128694e-07, 1.77574279496e-07,
  1.83813550477e-07, 1.92166040885e-07, 2.05295471952e-07, 2.22600839893e-07
};
