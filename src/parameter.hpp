#pragma once

#include <string>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <array>

#include "mymath.hpp"
#include "defs.hpp"

class Parameter {
  void LoadMacroParam() {
    const std::string fname = cur_dir + "/macro_IF.txt";
    std::ifstream fin(fname.c_str());
    CHECK_FILE_OPEN(fin);

    fin >> tempera >> sys_size >> headN >> tailN >> dt >> grid_leng.x >> L.x >> L.y >> L.z
	>> coef_prob[0] >> coef_prob[1] >> coef_prob[2] >> prob_cutof
	>> ls_grid_.x >> ls_grid_.y >> ls_grid_.z >> ls_lambda
	>> all_time_ >> time_step_mic_ >> time_step_mac_ >> time_step_vt_ >> chem_beg_;

    CHECK_FILESTREAM_IS_OK(fin);
    CHECK_FILE_IS_EOF(fin);
    
    grid_leng.z = grid_leng.y = grid_leng.x;
    
    ampN = (headN < tailN) ? headN : tailN;

    if (HYPHIL_N > REAC_PART) {
      hN = REAC_PART * headN + (HYPHIL_N - REAC_PART) * tailN;
      bN = HYPHOB_N * tailN;
    } else {
      hN = HYPHIL_N * headN;
      bN = (HYPHOB_N - (REAC_PART - HYPHIL_N)) * tailN + headN * (REAC_PART - HYPHIL_N);
    }

    wN = sys_size - bN - hN;

    iL	      = 1.0 / L;
    hL	      = 0.5 * L;
    ihL	      = 2.0 / L;

    dt_c      = dt * static_cast<double>(COL_FREQ);
    inv_dt_sq = 1.0 / std::sqrt(dt);
    
    for (int i = 0; i < 3; i++) {
      grid_numb[i]    = static_cast<int>(L[i] / grid_leng[i]);
      ls_grid_num_[i] = static_cast<int>(L[i] / ls_grid_[i]);
    }
      
    for (int i = 0; i < 3; i++) {
      grid_leng[i] = L[i] / grid_numb[i];
      ls_grid_[i]  = L[i] / ls_grid_num_[i];
    }
    all_grid = grid_numb[0] * grid_numb[1] * grid_numb[2];
    
    i_grid_leng = 1.0 / grid_leng;
    i_ls_grid_  = 1.0 / ls_grid_;
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
    
    CHECK_FILESTREAM_IS_OK(fin);
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
    CHECK_FILESTREAM_IS_OK(fin);
    CHECK_FILE_IS_EOF(fin);
  }

  void CheckMacroParam() const {
    CHECK_EQUATION(sys_size > 0, sys_size);
    
    CHECK_EQUATION(all_time_ > 0, all_time_);
    CHECK_EQUATION(time_step_mic_ > 0, time_step_mic_);
    CHECK_EQUATION(time_step_mac_ > 0, time_step_mac_);
    CHECK_EQUATION(time_step_vt_ > 0, time_step_vt_);
    CHECK_EQUATION(chem_beg_ > 0, chem_beg_);

    CHECK_EQUATION(wN > 0, wN);
    CHECK_EQUATION(hN > 0, hN);
    CHECK_EQUATION(bN > 0, bN);
    CHECK_EQUATION(ampN > 0, ampN);

    CHECK_EQUATION(tailN > 0, tailN);
    CHECK_EQUATION(headN > 0, headN);
    
    CHECK_EQUATION(std::isfinite(dt), dt);
    CHECK_EQUATION(std::isfinite(tempera), tempera);
    CHECK_EQUATION(std::isfinite(inv_dt_sq), inv_dt_sq);
    CHECK_EQUATION(std::isfinite(dt_c), dt_c);
    CHECK_EQUATION(std::isfinite(ls_lambda), ls_lambda);
    
    CHECK_EQUATION(L > double3(0.0), L);
    CHECK_EQUATION(iL > double3(0.0), iL);
    CHECK_EQUATION(hL > double3(0.0), hL);
    CHECK_EQUATION(ihL > double3(0.0), ihL);
    CHECK_EQUATION(grid_leng > double3(0.0), grid_leng);
    CHECK_EQUATION(i_grid_leng > double3(0.0), i_grid_leng);
    CHECK_EQUATION(ls_grid_ > double3(0.0), ls_grid_);
    CHECK_EQUATION(i_ls_grid_ > double3(0.0), i_ls_grid_);

    CHECK_EQUATION(ls_grid_num_[0] > 0, ls_grid_num_[0]);
    CHECK_EQUATION(ls_grid_num_[1] > 0, ls_grid_num_[1]);
    CHECK_EQUATION(ls_grid_num_[2] > 0, ls_grid_num_[2]);

    CHECK_EQUATION(grid_numb[0] > 0, grid_numb[0]);
    CHECK_EQUATION(grid_numb[1] > 0, grid_numb[1]);
    CHECK_EQUATION(grid_numb[2] > 0, grid_numb[2]);
    CHECK_EQUATION(all_grid > 0, all_grid);
    
    CHECK_EQUATION(ls_grid_[0] > 0.0, ls_grid_[0]);
    CHECK_EQUATION(ls_grid_[1] > 0.0, ls_grid_[1]);
    CHECK_EQUATION(ls_grid_[2] > 0.0, ls_grid_[2]);

    CHECK_EQUATION(i_ls_grid_[0] > 0.0, i_ls_grid_[0]);
    CHECK_EQUATION(i_ls_grid_[1] > 0.0, i_ls_grid_[1]);
    CHECK_EQUATION(i_ls_grid_[2] > 0.0, i_ls_grid_[2]);
    
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

  enum {
    HYPHIL_N = 1,
    HYPHOB_N = 3,
    ALL_UNIT_N = HYPHIL_N + HYPHOB_N,
    REAC_PART = 1,
    TAIL_PART = ALL_UNIT_N - REAC_PART,
    BUF_SIZE = 400,
    COL_FREQ = 10,
    EQUIL_TIME = 100000,
    LS_EQUIL_TIME = 5000,
  };

  static constexpr double b_leng   = 0.5;
  static constexpr double i_bleng  = 1.0 / b_leng;
  static constexpr double ch_leng  = 0.69;

  static int sys_size;
  
  int all_time_ = -1, time_step_mic_ = -1, time_step_mac_ = -1, time_step_vt_ = -1, chem_beg_ = -1;
  
  int wN = -1, bN = -1, hN = -1, ampN = -1, tailN = -1, headN = -1;
  double dt		= std::numeric_limits<double>::signaling_NaN();
  double tempera	= std::numeric_limits<double>::signaling_NaN();
  double inv_dt_sq	= std::numeric_limits<double>::signaling_NaN();
  double dt_c		= std::numeric_limits<double>::signaling_NaN();
  double ls_lambda      = std::numeric_limits<double>::signaling_NaN();
  
  double3 L, iL, hL, ihL, grid_leng, i_grid_leng;
  double3 ls_grid_, i_ls_grid_;
  
  int grid_numb[3] = {-1, -1, -1}, all_grid = -1;
  std::array<int, 3> ls_grid_num_;

  float coef_prob[3] = {-1.0, -1.0, -1.0}, prob_cutof = -1.0;

  explicit Parameter(char* cur_dir_) {
    cur_dir = cur_dir_;

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
	intrparams.cf_repul[i][j] = intrparams.cf_gamma[i][j] = intrparams.cf_sigma[i][j] = std::numeric_limits<double>::signaling_NaN();
    
    intrparams.cf_bend	 = std::numeric_limits<double>::signaling_NaN();
    intrparams.cf_spring = std::numeric_limits<double>::signaling_NaN();
    
    binfo.bind_center.x	= binfo.bind_center.y = binfo.bind_center.z = std::numeric_limits<double>::signaling_NaN();
    binfo.bind_radius	= std::numeric_limits<double>::signaling_NaN();
    binfo.bind_coef	= std::numeric_limits<double>::signaling_NaN();

    ls_grid_num_ = {-1, -1, -1};
  }

  Interactions GetIntractions() const {
    return intrparams;
  }
  
  BindInfo GetBindInform() const {
    return binfo;
  }

  int all_time() const {
    return all_time_;
  }
  
  int time_step_mic() const {
    return time_step_mic_;
  }
  
  int time_step_mac() const {
    return time_step_mac_;
  }

  int time_step_vt() const {
    return time_step_vt_;
  }

  int chem_beg() const {
    return chem_beg_;
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

    ost << "bind information" << std::endl;
    PRT_WITH_TAG(binfo.bind_center);
    PRT_WITH_TAG(binfo.bind_radius);
    PRT_WITH_TAG(binfo.bind_coef);

    ost << "interactions" << std::endl;
    for(int i = 0; i < 3; i++)
      ost << "sigma " << intrparams.cf_sigma[i][0] << " " << intrparams.cf_sigma[i][1] << " " << intrparams.cf_sigma[i][2] << std::endl;
    for(int i = 0; i < 3; i++)
      ost << "gamma " << intrparams.cf_gamma[i][0] << " " << intrparams.cf_gamma[i][1] << " " << intrparams.cf_gamma[i][2] << std::endl;
    for(int i = 0; i < 3; i++)
      ost << "repul " << intrparams.cf_repul[i][0] << " " << intrparams.cf_repul[i][1] << " " << intrparams.cf_repul[i][2] << std::endl;
    PRT_WITH_TAG(intrparams.cf_spring);
    PRT_WITH_TAG(intrparams.cf_bend);

    PRT_WITH_TAG(HYPHIL_N);
    PRT_WITH_TAG(HYPHOB_N);
    PRT_WITH_TAG(ALL_UNIT_N);
    PRT_WITH_TAG(REAC_PART);
    PRT_WITH_TAG(TAIL_PART);
    PRT_WITH_TAG(BUF_SIZE);
    PRT_WITH_TAG(COL_FREQ);
    PRT_WITH_TAG(EQUIL_TIME);
    
    PRT_WITH_TAG(b_leng);
    PRT_WITH_TAG(i_bleng);
    PRT_WITH_TAG(ch_leng);
    PRT_WITH_TAG(sys_size);

    const double rho = sys_size / (L.x * L.y * L.z);
    PRT_WITH_TAG(rho);
    
    PRT_WITH_TAG(wN); PRT_WITH_TAG(bN); PRT_WITH_TAG(hN); PRT_WITH_TAG(ampN);
    PRT_WITH_TAG(tailN); PRT_WITH_TAG(headN);

    PRT_WITH_TAG(dt); PRT_WITH_TAG(tempera); PRT_WITH_TAG(inv_dt_sq);
    PRT_WITH_TAG(dt_c); PRT_WITH_TAG(ls_lambda);
    
    PRT_WITH_TAG(L); PRT_WITH_TAG(iL); PRT_WITH_TAG(hL); PRT_WITH_TAG(ihL);
    PRT_WITH_TAG(grid_leng); PRT_WITH_TAG(i_grid_leng);
    PRT_WITH_TAG(ls_grid_); PRT_WITH_TAG(i_ls_grid_);
    
    PRT_WITH_TAG(grid_numb[0]); PRT_WITH_TAG(grid_numb[1]); PRT_WITH_TAG(grid_numb[2]);
    PRT_WITH_TAG(all_grid);
    PRT_WITH_TAG(ls_grid_num_[0]); PRT_WITH_TAG(ls_grid_num_[1]); PRT_WITH_TAG(ls_grid_num_[2]);
    
    PRT_WITH_TAG(coef_prob[0]);
    PRT_WITH_TAG(coef_prob[1]);
    PRT_WITH_TAG(coef_prob[2]);
    PRT_WITH_TAG(prob_cutof);

    ost << "COMPILE INFORM\n";
#ifdef NO_THERMO_STAT
    ost << "NO_THERMO_STAT is defined.\n";
#endif
#ifdef Y_REFLECT_BOUND
    ost << "Y_REFLECT_BOUND is defined.\n";
#endif
#ifdef ADD_BIND_POTENT
    ost << "ADD_BIND_POTENT is defined.\n";
#endif
#ifdef ADD_POSRES
    ost << "ADD_POSRES is defined.\n";
#endif
#ifdef CALC_HEIGHT
    ost << "CALC_HEIGHT is defined.\n";
#endif
#ifdef VISC_KERN_SQRT
    ost << "VISC_KERN_SQRT is defined.\n";
#endif
#ifdef CALC_LOC_STRESS
    ost << "CALC_LOC_STRESS is defined.\n";

#ifdef AT_FORCE_CENTER
    ost << "AT_FORCE_CENTER is defined.\n";
#endif
#ifdef AT_HYPOT
    ost << "AT_HYPOT is defined.\n";
#endif
#ifdef AT_MSP_GM
    ost << "AT_MSP_GM is defined.\n";
#endif
#ifdef AT_MSP_LM
    ost << "AT_MSP_LM is defined.\n";
#endif

#endif // end of CALC_LOC_STRESS

#undef PRT_WITH_TAG
  }

  void DumpGraphicDat(const int all_time, const int time_step) const {
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
