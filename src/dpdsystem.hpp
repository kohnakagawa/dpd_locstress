#pragma once

#include <fstream>
#include "parameter.hpp"

class Initializer;
class Observer;
class B_sorter;
class F_calculator;
class ChemManager;
class RNG;

#ifdef PROFILE_MODE
class Profile;
#endif

class dpdsystem {
  const Parameter       *p_param = nullptr;
  Observer		*p_observer = nullptr;

#ifdef PROFILE_MODE
  Profile               *p_profile = nullptr;
#endif

#if 0
#warning "memory alignment is 32 bit."
  //memory alignment
  enum {
    MEM_ALIGNMENT = 32,
  };
#endif

  //for bucket sort
  CellList             celllist;

  //for chemical reaction
  ChemInfo             cheminfo;

  void          AllocMem(const Parameter& param);
  void          DeallocMem();
  
  void          VelocVerlet1();
  void          VelocVerlet2();
  
  void          ReflectBound(double temp_r, double& r, double3& v, double3& F, double L) {
    if (temp_r >= L) {
      const double dw2p = temp_r - L;
      r   =  temp_r - 2.0 * dw2p;
      v  *= -1.0;
      F  *= -1.0;
    } else if (temp_r < 0.0) {
      r   = -temp_r;
      v  *= -1.0;
      F  *= -1.0;
    } else {
      r   = temp_r;
    }
  }

  void     RemoveCMDrift();
#ifdef ADD_POSRES
  void     SetPosresPrtcl();
#endif
  int      ReadParticleConfig();
  void     SetBondedParameter();
  void     CheckInitialized() const;
  int      GenParticles(ChemManager& chemmanage);
  
  void     FirstStep(B_sorter&, ChemManager&, F_calculator&);
public:
  //particle data
  double3	*pr = nullptr, *pv = nullptr, *delta_sumr = nullptr, *pv_bef = nullptr, *force = nullptr, *force_bef = nullptr;
  par_prop	*prop = nullptr;

#ifdef ADD_POSRES
  double3 *pr_base = nullptr;
  bool    *rest_on = nullptr;
#endif
  
  explicit dpdsystem(const Parameter& p_param_);
  ~dpdsystem();

  void Initialize();
  

  void     Execute(const int all_time,
		   const int time_step_mic,
		   const int time_step_mac,
		   const int time_step_vt,
		   const int chem_beg);

  int      GetPrtclIdx   (int i) const {   return celllist.prtcl_idx[i];    }
  int      GetLipidElemIdx(int i) const {  return cheminfo.lip_elem_idx[i];}
  bool     GetLipidChemConf(int i) const { return cheminfo.lipid_chem[i]; }
  int      GetPartIdx(int i) const { return cheminfo.part_idx[i]; }

  void     SetMoment  (int i, const double3& v ){ pv[i]	  = v; pv_bef[i] = v;}
  void     SetPosition(int i, const double3& r ){ pr[i]	  = r;   }
  void     SetProp    (int i, const par_prop& p){ prop[i]	  = p; }
  void     SetChemConf(int i, bool f)     { cheminfo.prtcl_chem[i] = f;  }
  void     SetLipidChem(int i, bool f)    { cheminfo.lipid_chem[i] = f;  }
  void     SetLipidIdx(int i, int idx)    { cheminfo.lipid_idx[i]  = idx;  }
  void     SetPartIdx(int i, int idx)     { cheminfo.part_idx[i]	  = idx;     }
  
  void     SetLipidUnit(int i, int unit)  { cheminfo.lipid_unit[i]    = unit; }
  void     SetLipidElemIdx(int i, int idx) {cheminfo.lip_elem_idx[i]  = idx;  }
  void     SetHeadElemIdx(int i, int idx){  cheminfo.head_elem_idx[i] = idx;  }
  void     SetTailElemIdx(int i, int idx){  cheminfo.tail_elem_idx[i] = idx;  }
};
