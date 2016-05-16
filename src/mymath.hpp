#pragma once
#include <cmath>
#include <cassert>

#define FAST_INV_SQRT

template <typename T>
T pow_n(T x, int n) {
  T a = 1.;
  
  while(n > 0){
    if(n % 2 == 0){
      x = x*x;
      n = n>>1;
    }else{
      a = a*x;
      n--;
    }
  }
  return a;
}

template<typename T>
T Pow_1_n(const T &a, const int n)
{
  const int  max = 1000;
  const T  eps = 1.0e-8;
  T  b = exp(log(a)/n);
  T revised_b;
  const T inv_n = 1.0 / n;
  for(int i=0; i<max; i++)
    {
      revised_b = b*(1.0 - inv_n)+a / (n*pow_n(b,n-1));
      if(fabs(b-revised_b)<eps) break;
      b = revised_b;
    }
  return revised_b;
}

//---Algorithm double(IEEE754)用------
template<typename T>
T rsqrtD(const T& x) 
{
#ifdef DEBUG
  assert(x > 0.);  
#endif
  const double         xHalf = 0.5 * (double)x;
  const long long int  tmp   = 0x5FE6EB50C7B537AAl - ( *(long long int*)&x >> 1);//initial guess
  double         xRes  = * (double*)&tmp;

  xRes *= ( 1.5 - ( xHalf * xRes * xRes ) );
#ifndef FAST_INV_SQRT
  xRes *= ( 1.5 - ( xHalf * xRes * xRes ) );
#endif
  return (T)xRes;
}

#undef FAST_INV_SQRT
