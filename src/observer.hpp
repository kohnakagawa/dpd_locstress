#pragma once
#include "parameter.hpp"
#include "dpdsystem.hpp"
#include <memory>

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
    LOC_STRESS,
    DECOMP_ERROR,
    
    NUM_FILE,
  };
  
  FILE* fp[NUM_FILE] = {nullptr};

  std::string Type2Fname(const int type, const Parameter& param);

  std::vector<double> loc_tempera, loc_dense;
  std::vector<double3> loc_vel, tail_cm_pos;
  std::unique_ptr<tensor3d[]> loc_stress_sum;
  int cnt_ls = 0;
  
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
  
  void DumpPrtclConfig(const dpdsystem &sDPD, const ChemInfo& cheminfo, const int time, FILE* fp);
  
public:
  explicit Observer(const Parameter& param);
  ~Observer();
  void Initialize(const Parameter& param);
  void DumpMacroVal(const dpdsystem& sDPD, const Parameter& param);
  void DumpPressure(const dpdsystem& sDPD, const Parameter& param, const tensor3d& vil);
  void DumpConfigTempera(const double configT);
  void DumpLocalVal(const dpdsystem &sDPD, const Parameter& param);
  void DumpTranject(const dpdsystem& sDPD, const ChemInfo& cheminfo, const int time);
  void DumpFinalConfig(const dpdsystem &sDPD, const ChemInfo& cheminfo, const int time);
  void DumpMembHeight(const dpdsystem& sDPD, const Parameter& param, const int time);
  void DumpForceDecompError(const double err);
  void AddLocalStress(const dpdsystem& sDPD, const Parameter& param, const std::vector<std::vector<tensor3d> >& buf_lstress);
  void DumpLocalStress(const Parameter& param);
};
