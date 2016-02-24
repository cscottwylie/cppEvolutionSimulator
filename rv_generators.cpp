#include "rv_generators.hpp"
#include "assert.h"
#include <cmath>
#include <iostream>

namespace evolve{

int rnd_binomial(double pp, int xn) {       // pp =probability heads, xn= # flips
  int j,n;
  double am,em,g,angle,p,bnl,sq,t,y;
  static double xnold=(-1.0),pold=(-1.0),pc,plog,pclog,oldg;
  p=(pp <= 0.5 ? pp : 1.0-pp);
  am=xn*p;
  if (xn < 25.0) {
    n=((int)(2.0*(xn)) + 1)/2;
    bnl=0.0;
    for (j=1;j<=n;j++)
      if (drand48() < p) bnl += 1.0;
  } else if (am < 1.0) {
    n=((int)(2.0*(xn)) + 1)/2;
    g=exp(-am);
    t=1.0;
    for (j=0;j<=n;j++) {
      t *= drand48();
      if (t < g) break;
    }
    bnl=(j <= n ? j : n);
  } else {
    if (xn != xnold) {
      oldg=lgamma(xn+1.0);
      xnold=xn;
    } if (p != pold) {
      pc=1.0-p;
      plog=log(p);
      pclog=log(pc);
      pold=p;
    }
    sq=sqrt(2.0*am*pc);
    do {
      do {
	angle=3.14159265358979323846*drand48();
	y=tan(angle);
	em=sq*y+am;
      } while (em < 0.0 || em >= (xn+1.0));
      em=floor(em);
      t=1.2*sq*(1.0+y*y)*exp(oldg-lgamma(em+1.0)
			     -lgamma(xn-em+1.0)+em*plog+(xn-em)*pclog);
    } while (drand48() > t);
    bnl=em;
  }
  if (p != pp) bnl=xn-bnl;
  return (int) bnl;
};

double rnd_gaussian(double mean, double stdev) {

  static int iset=0;
  static double gset;

  double fac, rsq, v1, v2;

  if( iset== 0) {
    do {
      v1= 2.0* drand48()- 1.0;
      v2= 2.0* drand48()- 1.0;

      rsq= v1*v1 + v2*v2;

    } while (rsq>= 1.0 || rsq== 0.0);

    fac= sqrt( -2.0* log( rsq) / rsq);
    gset= v1* fac;
    iset= 1;

    return mean+ stdev* v2* fac;

  } else {
    iset= 0;
    return mean+ stdev* gset;
  }
}

double rnd_konstantine(){                // analytical dG distribution from PNAS '07
                       // probability    // deltaG
  const double weights[12]= { 0.0000,   // -20.0000
    					      0.0031,   // -18.1818
    					      0.0093,   // -16.3636
    					      0.0204,   // -14.5455
   					          0.0388,   // -12.7273
    					      0.0665,   // -10.9091
    				    	  0.1047,   // -9.0909
    					      0.1516,   // -7.2727
    					      0.1984,   // -5.4545
    					      0.2236,   // -3.6364
    					      0.1836,   // -1.8182
         					  0         // 0            /// these probabilities should sum to 1
												};
  const double dG[12]=      {-20.0000,                  // evenly spaced dG b/w -20, 0
   						     -18.1818,
  	  					     -16.3636,
  						     -14.5455,
 						     -12.7273,
  						     -10.9091,
  						     -9.0909,
  	   					     -7.2727,
   						     -5.4545,
   						     -3.6364,
   						     -1.8182,
       					     0      };

  double ch= rnd_uniform();
  int i= 11;
  while( ( ch+= weights[ i] ) < 1) {
    --i;
  };
  assert( i>= 0);
  assert( i<= 11);

  return dG[ i];
  }


}


