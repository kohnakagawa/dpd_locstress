#pragma once
 
#include <fstream>
#include <string>
#include <sys/time.h>
 
class Profile{
  std::ofstream fout;
  timeval beg;
  double dur_t[4];
  
  inline double time_dif(timeval &beg,timeval &end){
    double delta_t;
    delta_t = end.tv_sec - beg.tv_sec + (end.tv_usec - beg.tv_usec) * 1.0e-6;
    return delta_t;
  }
public:
  enum{
    B_SORT=0,
    TIME_EVOLVE=1,
    F_CALC=2,
    DET_CHEM=3
  };

  explicit Profile(const char* cur_dir){
    for(int i=0; i<4; i++) dur_t[i] = 0.0;
    std::string temp_st = cur_dir;
    temp_st += "/elapse_time.txt";
    fout.open(temp_st.c_str());
  }
  ~Profile(){
#define PRT(arg1,arg2) fout << arg1 << " " << arg2 << std::endl;
    PRT("sort       ",dur_t[B_SORT]);
    PRT("time evolve",dur_t[TIME_EVOLVE]);
    PRT("force calc ",dur_t[F_CALC]);
    PRT("chem part  ",dur_t[DET_CHEM]);
#undef PRT
  }
  inline void w_start(){
    gettimeofday(&beg,NULL);
  }
  inline void w_stop(int kind){
    timeval end;
    gettimeofday(&end,NULL);
    dur_t[kind] += time_dif(beg,end);
  }
};
