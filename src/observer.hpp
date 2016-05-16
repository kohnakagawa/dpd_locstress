#pragma once
#include "parameter.hpp"
#include "dpdsystem.hpp"

class Observer {
  // observer type
  enum {
    MACRO_VAL = 0,
    PTCL_CONFIG,
    LOCAL_VAL,
    PRESSURE,
    GYRATION,
    CONFIG_TEMP,
    FINAL_CONFIG,
    HEIGHT_DIST,
    
    NUM_FILE,
  };
  
  FILE* fp[NUM_FILE] = {nullptr};

  std::string Type2Fname(const int type, const Parameter& param);

  std::vector<double> loc_tempera, loc_dense;
  std::vector<double3> loc_vel, tail_cm_pos;
  
  void MinImage(double3& a, const Parameter& param) {
    a.x -= param.L.x * std::round(a.x * param.iL.x);
    a.y -= param.L.y * std::round(a.y * param.iL.y);
    a.z -= param.L.z * std::round(a.z * param.iL.z);
  }

  int GenHash(const double3 &r,const Parameter& param){
    int i = static_cast<int>(r.x * param.i_grid_leng.x);
    int j = static_cast<int>(r.y * param.i_grid_leng.y);
    int k = static_cast<int>(r.z * param.i_grid_leng.z);

    if (i == param.grid_numb[0]) i--;
    if (j == param.grid_numb[1]) j--;
    if (k == param.grid_numb[2]) k--;
  
    const int hash = i + param.grid_numb[0] * (j + k * param.grid_numb[1]);
#ifdef DEBUG
    assert(0 <= hash);
    assert(hash < param.all_grid);
#endif
    return hash;
  }

  //ASSUME: membrane is square shape.
  double cut_r = 2.0, i_cut_r = -1.0;
  int cut_grid = -1;
  std::vector<double> hei; std::vector<int> hei_elem;
  
  int GenHash2d(const double3& r) {
    int ix = static_cast<int>(r.x * i_cut_r);
    int iz = static_cast<int>(r.z * i_cut_r);
    if (ix == cut_grid) ix--;
    if (iz == cut_grid) iz--;
    return iz + cut_grid * ix;
  }
  
  double	CalcKinTempera(const dpdsystem& sDPD);
  double	CalcDiffs(const dpdsystem& sDPD);
  double	CalcGyration(const dpdsystem& sDPD,const Parameter& param);
  double	CalcThickness(const dpdsystem& sDPD,const Parameter& param);
  double	CalcOrientOrder(const dpdsystem& sDPD,const Parameter& param);  
  void          CalcMembraneHeight(const dpdsystem& sDPD);
  
  void DumpPrtclConfig(const dpdsystem &sDPD, const ChemInfo& cheminfo, FILE* fp);
  
public:
  explicit Observer(const Parameter& param);
  ~Observer();
  void Initialize(const Parameter& param);
  void DumpMacroVal(const dpdsystem& sDPD, const Parameter& param);
  void DumpPressure(const dpdsystem& sDPD, const Parameter& param, const double3& vil);
  void DumpConfigTempera(const double configT);
  void DumpLocalVal(const dpdsystem &sDPD, const Parameter& param);
  void DumpTranject(const dpdsystem& sDPD, const ChemInfo& cheminfo);
  void DumpFinalConfig(const dpdsystem &sDPD, const ChemInfo& cheminfo);
  void DumpMembHeight(const dpdsystem& sDPD, const Parameter& param, const int time);
};
