#pragma once

#include <string>
#include <fstream>
#include <cstdlib>
#include <limits>

#include "mymath.hpp"
#include "defs.hpp"

class Parameter {
  void LoadMacroParam() {
    const std::string fname = cur_dir + "/macro_IF.txt";
    std::ifstream fin(fname.c_str());
    CHECK_FILE_OPEN(fin);

    double rho = 0.0;
    fin >> rho >> tempera >> headN >> tailN >> dt >> grid_leng.x >> L.y >> coef_prob[0] >> coef_prob[1] >> coef_prob[2] >> prob_cutof;
    L.z = L.x = std::sqrt(SYS_SIZE / (rho * L.y));
    grid_leng.z = grid_leng.y = grid_leng.x; //ASSUME: grid_leng.x = grid_leng.y = grid_leng.z
    
    ampN = (headN < tailN) ? headN : tailN;
    
    CHECK_FILE_IS_EOF(fin);

    if (HYPHIL_N > REAC_PART) {
      hN = REAC_PART * headN + (HYPHIL_N - REAC_PART) * tailN;
      bN = HYPHOB_N * tailN;
    } else {
      hN = HYPHIL_N * headN;
      bN = (HYPHOB_N - (REAC_PART - HYPHIL_N)) * tailN + headN * (REAC_PART - HYPHIL_N);
    }

    wN = SYS_SIZE - bN - hN;

    iL	      = 1.0 / L;
    hL	      = 0.5 * L;
    ihL	      = 2.0 / L;

    dt_c      = dt * static_cast<double>(COL_FREQ);
    inv_dt_sq = 1.0 / std::sqrt(dt);
    
    for (int i = 0; i < 3; i++)
      grid_numb[i] = static_cast<int>(L[i] / grid_leng[i]);
    for (int i = 0; i < 3; i++)
      grid_leng[i] = L[i] / grid_numb[i];
    
    all_grid = grid_numb[0] * grid_numb[1] * grid_numb[2];
    i_grid_leng = 1.0 / grid_leng;
  }
  
  void CalculateGammaWithHarmonicMean(const int i, const int j) {
    if (intrparams.cf_sigma[i][j] < 0.0) {
      intrparams.cf_gamma[i][j] = intrparams.cf_gamma[j][i] = 2.0 / ( (1.0 / intrparams.cf_gamma[i][i]) + (1.0 / intrparams.cf_gamma[j][j]) );
      intrparams.cf_sigma[i][j] = intrparams.cf_sigma[j][i] = std::sqrt(2.0 * intrparams.cf_gamma[i][j] * tempera);
    }
  }
  
  void LoadMicroParam() {
    std::string fname = cur_dir + "/micro_IF.txt";
    std::ifstream fin(fname.c_str());
    CHECK_FILE_OPEN(fin);

    fin >> intrparams.cf_sigma[0][0] >> intrparams.cf_sigma[0][1] >> intrparams.cf_sigma[0][2] 
	>> intrparams.cf_sigma[1][1] >> intrparams.cf_sigma[1][2] >> intrparams.cf_sigma[2][2] 
	>> intrparams.cf_spring      >> intrparams.cf_bend
	>> intrparams.cf_repul[0][0] >> intrparams.cf_repul[0][1] >> intrparams.cf_repul[0][2]
	>> intrparams.cf_repul[1][1] >> intrparams.cf_repul[1][2] >> intrparams.cf_repul[2][2];
    
    CHECK_FILE_OPEN(fin);

    intrparams.cf_repul[1][0] = intrparams.cf_repul[0][1];
    intrparams.cf_repul[2][0] = intrparams.cf_repul[0][2];
    intrparams.cf_repul[2][1] = intrparams.cf_repul[1][2];

    intrparams.cf_sigma[1][0] = intrparams.cf_sigma[0][1];
    intrparams.cf_sigma[2][0] = intrparams.cf_sigma[0][2];
    intrparams.cf_sigma[2][1] = intrparams.cf_sigma[1][2];

    intrparams.cf_gamma[0][0] = 0.5 * intrparams.cf_sigma[0][0] * intrparams.cf_sigma[0][0] / tempera;
    intrparams.cf_gamma[1][1] = 0.5 * intrparams.cf_sigma[1][1] * intrparams.cf_sigma[1][1] / tempera;
    intrparams.cf_gamma[2][2] = 0.5 * intrparams.cf_sigma[2][2] * intrparams.cf_sigma[2][2] / tempera;
    intrparams.cf_gamma[0][1] = 0.5 * intrparams.cf_sigma[0][1] * intrparams.cf_sigma[0][1] / tempera;
    intrparams.cf_gamma[0][2] = 0.5 * intrparams.cf_sigma[0][2] * intrparams.cf_sigma[0][2] / tempera;
    intrparams.cf_gamma[1][2] = 0.5 * intrparams.cf_sigma[1][2] * intrparams.cf_sigma[1][2] / tempera;

    intrparams.cf_gamma[1][0] = intrparams.cf_gamma[0][1];
    intrparams.cf_gamma[2][0] = intrparams.cf_gamma[0][2];
    intrparams.cf_gamma[2][1] = intrparams.cf_gamma[1][2];

    //NOTE: Water-head interaction is determined by harmonic mean.
    //      This is employed when intrparams.cf_sigma[i][j] < 0.0.
    CalculateGammaWithHarmonicMean(0, 1);
    CalculateGammaWithHarmonicMean(0, 2);
    CalculateGammaWithHarmonicMean(1, 2);
  }
  
  void LoadBindParam() {
    std::string fname = cur_dir + "/bind_info.txt";
    std::ifstream fin(fname.c_str());
    CHECK_FILE_OPEN(fin);
    fin >> binfo.bind_radius >> binfo.bind_coef >> binfo.bind_center.x >> binfo.bind_center.y >> binfo.bind_center.z;
    CHECK_FILE_IS_EOF(fin);
  }

  void CheckMacroParam() const {
    CHECK_EQUATION(wN > 0, wN);
    CHECK_EQUATION(hN > 0, hN);
    CHECK_EQUATION(bN > 0, bN);
    CHECK_EQUATION(ampN > 0, ampN);

    CHECK_EQUATION(tailN > 0, tailN);
    CHECK_EQUATION(headN > 0, headN);
    
    CHECK_EQUATION(std::isfinite(dt), dt);
    CHECK_EQUATION(std::isfinite(tempera), tempera);
    CHECK_EQUATION(std::isfinite(inv_dt_sq), inv_dt_sq);
    
    CHECK_EQUATION(L.isfinite3(), L);
    CHECK_EQUATION(iL.isfinite3(), iL);
    CHECK_EQUATION(hL.isfinite3(), hL);
    CHECK_EQUATION(grid_leng.isfinite3(), grid_leng);
    CHECK_EQUATION(i_grid_leng.isfinite3(), i_grid_leng);

    CHECK_EQUATION(std::isfinite(grid_numb[0]), grid_numb[0]);
    CHECK_EQUATION(std::isfinite(grid_numb[1]), grid_numb[1]);
    CHECK_EQUATION(std::isfinite(grid_numb[2]), grid_numb[2]);
    CHECK_EQUATION(std::isfinite(all_grid), all_grid);

    CHECK_EQUATION(coef_prob[0] > 0.0, coef_prob[0]);
    CHECK_EQUATION(coef_prob[1] > 0.0, coef_prob[1]);
    CHECK_EQUATION(coef_prob[1] < 1.0, coef_prob[1]);
    CHECK_EQUATION(coef_prob[2] > 0.0, coef_prob[2]);
    CHECK_EQUATION(prob_cutof   > 0.0, prob_cutof);
  }

  void CheckMicroParam() const {
    for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) {
	CHECK_EQUATION(intrparams.cf_sigma[i][j] >= 0.0, intrparams.cf_sigma[i][j]);
	CHECK_EQUATION(std::isfinite(intrparams.cf_sigma[i][j]), intrparams.cf_sigma[i][j]);
	CHECK_EQUATION(intrparams.cf_gamma[i][j] >= 0.0, intrparams.cf_gamma[i][j]);
	CHECK_EQUATION(std::isfinite(intrparams.cf_gamma[i][j]), intrparams.cf_gamma[i][j]);
	CHECK_EQUATION(intrparams.cf_repul[i][j] >= 0.0, intrparams.cf_repul[i][j]);
	CHECK_EQUATION(std::isfinite(intrparams.cf_repul[i][j]), intrparams.cf_repul[i][j]);
	
	CHECK_EQUATION(intrparams.cf_sigma[i][j] == intrparams.cf_sigma[j][i], intrparams.cf_sigma[i][j]);
	CHECK_EQUATION(intrparams.cf_gamma[i][j] == intrparams.cf_gamma[j][i], intrparams.cf_gamma[i][j]);
	CHECK_EQUATION(intrparams.cf_repul[i][j] == intrparams.cf_repul[j][i], intrparams.cf_repul[i][j]);
      }
    CHECK_EQUATION(intrparams.cf_bend >= 0.0, intrparams.cf_bend);
    CHECK_EQUATION(std::isfinite(intrparams.cf_bend), intrparams.cf_bend);
    CHECK_EQUATION(intrparams.cf_spring >= 0.0, intrparams.cf_spring);
    CHECK_EQUATION(std::isfinite(intrparams.cf_spring), intrparams.cf_spring);
  }

  void CheckBindInfo() const {
    CHECK_EQUATION(binfo.bind_center.x >= 0.0, binfo.bind_center.x);
    CHECK_EQUATION(binfo.bind_center.y >= 0.0, binfo.bind_center.y);
    CHECK_EQUATION(binfo.bind_center.z >= 0.0, binfo.bind_center.z);
    CHECK_EQUATION(binfo.bind_coef >= 0.0, binfo.bind_coef);
    CHECK_EQUATION(binfo.bind_radius >= 0.0, binfo.bind_radius);
  }

public:
  std::string cur_dir;
  BindInfo binfo;
  Interactions intrparams;

  //change system size
  enum {
    HYPHIL_N = 1,
    HYPHOB_N = 3,
    ALL_UNIT_N = HYPHIL_N + HYPHOB_N,
    REAC_PART = 1, //bond connect/divide 0~1
    TAIL_PART = ALL_UNIT_N - REAC_PART,
    SYS_SIZE = 331776,
    //normal box 48*48*48*3 = 331776,  40*40*40*3 = 192000, 36*36*36*3 = 139968, 32*32*32*3 = 98304, 28*28*28*3 = 65856, 24*24*24*3 = 41472, 20*20*20*3 = 24000 16*16*16*3=12288 12*12*12*3=5184 8*8*8*3=1536
    //cuboid box 32*36*32*3 = 110592, 32*40*32*3 = 122880, 32*44*32*3 = 135168 12*16*12*3 = 6912 36*36*48*3 = 186624
    BUF_SIZE = 400,
    COL_FREQ = 10,
    EQUIL_TIME = 500000,
  };

  static constexpr double lip_area = 2.0 * 36.0 * 36.0 / 4950.0; //NOTE: This value depends on molecular structure.
  static constexpr double b_leng   = 0.8;
  static constexpr double i_bleng  = 1.0 / b_leng; // i_bleng = 1.0 / b_leng 

  static constexpr double ch_leng  = 0.69;
  
  int wN = -1, bN = -1, hN = -1, ampN = -1, tailN = -1, headN = -1;
  double dt		= std::numeric_limits<double>::signaling_NaN();
  double tempera	= std::numeric_limits<double>::signaling_NaN();
  double inv_dt_sq	= std::numeric_limits<double>::signaling_NaN();
  double dt_c		= std::numeric_limits<double>::signaling_NaN();
  
  double3 L, iL, hL, ihL, grid_leng, i_grid_leng;

  int grid_numb[3] = {-1, -1, -1}, all_grid = -1;
  float coef_prob[3] = {-1.0, -1.0, -1.0}, prob_cutof = -1.0;

  explicit Parameter(char* cur_dir_) {
    cur_dir = cur_dir_;

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	intrparams.cf_repul[i][j] = intrparams.cf_gamma[i][j] = intrparams.cf_sigma[i][j] = std::numeric_limits<double>::signaling_NaN();
    
    intrparams.cf_bend		= std::numeric_limits<double>::signaling_NaN();
    intrparams.cf_spring	= std::numeric_limits<double>::signaling_NaN();
    
    binfo.bind_center.x		= binfo.bind_center.y = binfo.bind_center.z = std::numeric_limits<double>::signaling_NaN();
    binfo.bind_radius		= std::numeric_limits<double>::signaling_NaN();
    binfo.bind_coef		= std::numeric_limits<double>::signaling_NaN();
  }

  Interactions GetIntractions() const {
    return intrparams;
  }
  
  BindInfo GetBindInform() const {
    return binfo;
  }
  
  void LoadParam() {
    LoadMacroParam();
    LoadBindParam();
    LoadMicroParam();
  }

  void LoadCheck() {
    CheckMacroParam();
    CheckMicroParam();
    CheckBindInfo();
  }

  void DumpAllParam(std::ostream& ost) const {
#define PRT_WITH_TAG(arg) ost << #arg << " = " << arg << std::endl;
    PRT_WITH_TAG(cur_dir);
    
    PRT_WITH_TAG(ALL_UNIT_N);
    PRT_WITH_TAG(REAC_PART);
    PRT_WITH_TAG(TAIL_PART);
    
    PRT_WITH_TAG(wN); PRT_WITH_TAG(bN); PRT_WITH_TAG(hN); PRT_WITH_TAG(ampN);

    PRT_WITH_TAG(tailN); PRT_WITH_TAG(headN);

    PRT_WITH_TAG(dt); PRT_WITH_TAG(tempera); PRT_WITH_TAG(inv_dt_sq);
    PRT_WITH_TAG(L); PRT_WITH_TAG(iL); PRT_WITH_TAG(hL); PRT_WITH_TAG(ihL);
    PRT_WITH_TAG(grid_leng); PRT_WITH_TAG(i_grid_leng);
    PRT_WITH_TAG(grid_numb[0]); PRT_WITH_TAG(grid_numb[1]); PRT_WITH_TAG(grid_numb[2]);
    PRT_WITH_TAG(all_grid);
    
    ost << "chemical parameters " << std::endl;
    PRT_WITH_TAG(coef_prob[0]);
    PRT_WITH_TAG(coef_prob[1]);
    PRT_WITH_TAG(coef_prob[2]);
    PRT_WITH_TAG(prob_cutof);

    ost << "interactions" << std::endl;
    for(int i = 0; i < 3; i++)
      ost << "sigma " << intrparams.cf_sigma[i][0] << " " << intrparams.cf_sigma[i][1] << " " << intrparams.cf_sigma[i][2] << std::endl;
    for(int i = 0; i < 3; i++)
      ost << "gamma " << intrparams.cf_gamma[i][0] << " " << intrparams.cf_gamma[i][1] << " " << intrparams.cf_gamma[i][2] << std::endl;
    for(int i = 0; i < 3; i++)
      ost << "repul " << intrparams.cf_repul[i][0] << " " << intrparams.cf_repul[i][1] << " " << intrparams.cf_repul[i][2] << std::endl;

    PRT_WITH_TAG(intrparams.cf_spring);
    PRT_WITH_TAG(intrparams.cf_bend);

    PRT_WITH_TAG(b_leng);
    PRT_WITH_TAG(ch_leng);

    ost << "bind information" << std::endl;
    PRT_WITH_TAG(binfo.bind_center);
    PRT_WITH_TAG(binfo.bind_radius);
    PRT_WITH_TAG(binfo.bind_coef);

#undef PRT_WITH_TAG
  }

  void DumpGraphicDat(int all_time, int time_step) const {
    std::string str;
    str = cur_dir + "/macro_data.txt";
    std::ofstream fout(str.c_str());
#define PRT(arg) fout << arg << std::endl
    PRT(wN);
    PRT(hN + bN);
    PRT(L.x);
    PRT(0.5);
    PRT(all_time);
    PRT(time_step);
#undef PRT
  }
};

//unit test
#if 0
int main(int argc, char* argv[]) {
  Parameter param(argv[1],argv[2]);
  param.LoadParam();
  param.LoadCheck(std::cout);
  param.DumpGraphicDat(100,100);
}
#endif
