// This class defines some random number generation functions
// Most of the underlying code was borrowed from the BOOST library 
#ifndef __QOL_RANDOM__
#define __QOL_RANDOM__

#include <cmath>

#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64) || defined(__WIN32) || defined(__WIN64) && !defined(__CYGWIN__)
#include <boost/random.hpp>

static boost::lagged_fibonacci1279 cpmInternalRandGen;  

class Random : public boost::lagged_fibonacci1279 {
public:
  Random () {}

  ~Random () {}

  inline void seed (int seed) 
  { cpmInternalRandGen.seed (seed); }

  inline double uniform(double min=0.0,double max=1.0)
  { return cpmInternalRandGen()*(max-min)+min; }

  /// returns integer from min to max-1. Default returns 0 or 1.
  inline int uniform_int(int min=0,int max=2)
  { return (int)floor(cpmInternalRandGen()*(max-min)+min); }
};

#else

#include <stdio.h>
#include <stdlib.h>

class Random {
public:
  Random () {}
  ~Random() {}
  inline void seed (int seed) { srand (seed); }

  inline double uniform(double min=0.0,double max=1.0)
  { return rand()*(max-min)+min; }

  /// returns integer from min to max-1. Default returns 0 or 1.
  inline int uniform_int(int min=0,int max=2) 
  { return (int)floor(double (rand()*(max-min)+min)); }
};
#endif

#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/

