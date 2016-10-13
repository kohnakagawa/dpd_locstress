#pragma once

#define INLINE __attribute__((always_inline))

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <ctime>

#include <random>

class RNG {
  const size_t n_rng;
  std::mt19937* gen = nullptr;
  std::uniform_real_distribution<> uni_dist; // NOTE: range is [0,1].
  std::normal_distribution<> *normal_dist = nullptr; // NOTE: mean=0 sd=1.0

  // NOTE:
  // uniform_real_distribution does not have state.
  // (uniform_real_distribution::reset does nothing!)
  // However, normal_distribution have state variable as class member.
  // So, we prepare a container for each thread.

public:
  RNG(const size_t seed, const size_t th_numb) : n_rng(th_numb) {
    gen = new std::mt19937 [th_numb];
    for (size_t i = 0; i < th_numb; i++) gen[i].seed(seed + i * 10);
    normal_dist = new std::normal_distribution<> [th_numb];
  }

  explicit RNG(const size_t th_numb) : n_rng(th_numb) {
    gen = new std::mt19937 [th_numb];
    const auto base_seed = static_cast<size_t>(time(nullptr));
    for (size_t i = 0; i < th_numb; i++) gen[i].seed(base_seed + i * 10);
    normal_dist = new std::normal_distribution<> [th_numb];
  }

  ~RNG(){
    delete [] gen;
    delete [] normal_dist;
  }

  INLINE double Uniform(const int tid) {
    return uni_dist(gen[tid]);
  }
  INLINE double Uniform(const int tid, const double up, const double dw) {
    return dw + (up - dw) * uni_dist(gen[tid]);
  }

  INLINE double Normal(const int tid) {
    return normal_dist[tid](gen[tid]);
  }
  INLINE double Normal(const int tid, const double mean, const double sd) {
    return mean + normal_dist[tid](gen[tid]) * std::sqrt(sd);
  }
};

#undef INLINE

#if 0
#include <iostream>
#include <iomanip>
#include <omp.h>

int main(){
  RNG rng(100, 4);

  std::cout << std::setprecision(16);

  // for(int i=0; i<100000; i++) std::cout << rng.Uniform(0) << std::endl;
  //for(int i=0; i<100000; i++) std::cout << rng.Uniform(0,2.0,4.0) << std::endl;
  //for(int i=0; i<100000; i++) std::cout << rng.Normal(0)  << std::endl;
  //for(int i=0; i<100000; i++) std::cout << rng.Normal(0,1.0,4.0) << std::endl;
#pragma omp parallel for
  for (int i = 0; i < 1000000; i++) {
    const int tid = omp_get_thread_num();
    const auto ran = rng.Normal(tid);
    // const auto ran = rng.Uniform(tid);
#pragma omp critical
    {
      std::cout << ran << std::endl;
    }
  }
}
#endif


