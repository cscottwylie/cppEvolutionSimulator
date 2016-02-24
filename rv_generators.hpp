#ifndef _RV_GENERATORS_
#define _RV_GENERATORS_

#include <cmath>
#include <iostream>
#include <stdlib.h>

namespace evolve{

inline double rnd_uniform();
inline int    rnd_int( int N);
inline double rnd_expo( double lambda);
int           rnd_binomial( double pp, int xn);
double        rnd_gaussian( double mean, double stdev);
double        rnd_konstantine();  // deltaG drawn from PNAS '07 equilibrium distribution

inline double rnd_uniform()           {return drand48();                   };
//inline double rnd_uniform()           {return gsl_rng_uniform(BaseRand); };
inline int    rnd_int(int N)          {return (int)(drand48()*N);          };
inline double rnd_expo(double lambda) {return -log(rnd_uniform()) / lambda;};

}

#endif
