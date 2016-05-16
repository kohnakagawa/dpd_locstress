#include "dpdsystem.hpp"
#include "chem_manager.hpp"
#include "allocate.hpp"
#include "observer.hpp"
#include "initializer.hpp"
#include "bucket_sorter.hpp"
#include "force_calculator.hpp"
#include "rng.hpp"
#include "profile.hpp"

#include <string>
#include <omp.h>
#include <ctime>
#include <limits>
#include <numeric>

dpdsystem::dpdsystem(const Parameter& p_param_) {
  p_param = &p_param_;
}

dpdsystem::~dpdsystem() {
  DeallocMem();
}

void dpdsystem::Initialize() {
  AllocMem(*p_param);
  
  // set values
#pragma omp parallel for // for NUMA
  for (int i = 0; i < Parameter::SYS_SIZE; i++) {
    pr[i]			= double3(std::numeric_limits<double>::signaling_NaN());
    pv[i]			= double3(std::numeric_limits<double>::signaling_NaN());
    prop[i]			= Water;
    pv_bef[i]			= double3(std::numeric_limits<double>::signaling_NaN());
    force[i]			= double3(std::numeric_limits<double>::signaling_NaN());
    force_bef[i]		= double3(std::numeric_limits<double>::signaling_NaN());
    delta_sumr[i]		= double3(std::numeric_limits<double>::signaling_NaN());
    celllist.next_dest[i]	= -1;

#ifdef ADD_POSRES
    rest_on[i]			= false;
    pr_base[i]			= double3(std::numeric_limits<double>::signaling_NaN());
#endif
  }

  for (int i = 0; i < p_param->all_grid; i++) {
    celllist.buck_elem[i]			= -1;
    celllist.buck_addrs[i]			= -1;
  }
  celllist.buck_addrs[p_param->all_grid]	= 0;

  const int int_max = std::numeric_limits<int>::max();
#pragma omp parallel for // for NUMA
  for (int i = 0; i < Parameter::SYS_SIZE; i++) {
    cheminfo.lipid_idx[i]	= int_max;
    cheminfo.part_idx[i]	= int_max;
    cheminfo.lipid_unit[i]	= int_max;
    cheminfo.near_water[i]	= 0;
    cheminfo.near_info[i].dist	= 4.0;
    cheminfo.near_info[i].idx	= -1;
  }
  
  for (int i = 0; i < p_param->headN * Parameter::REAC_PART; i++)
    cheminfo.head_elem_idx[i] = int_max;

  for (int i = 0; i < p_param->tailN * Parameter::TAIL_PART; i++)
    cheminfo.tail_elem_idx[i] = int_max;
  
  for (int i = 0; i < p_param->ampN * Parameter::ALL_UNIT_N; i++)
    cheminfo.lip_elem_idx[i] = int_max;

  for (int i = 0; i < Parameter::SYS_SIZE; i++) {
    cheminfo.near_water[i]	= 0;
    cheminfo.near_info[i].dist	= 4.0;
    cheminfo.near_info[i].idx	= -1;
  }
}

void dpdsystem::AllocMem(const Parameter& param) {
  p_observer	= new Observer(param);

#ifdef PROFILE_MODE
  p_profile     = new Profile(p_param_.cur_dir.c_str());
#endif

  celllist.buck_elem   = new int      [param.all_grid    ];
  celllist.buck_addrs  = new int      [param.all_grid + 1];
  celllist.next_dest   = new int      [Parameter::SYS_SIZE];

#if 0
  pr         = allocate<double3,MEM_ALIGNMENT >(Parameter::SYS_SIZE);
  pv         = allocate<double3,MEM_ALIGNMENT >(Parameter::SYS_SIZE);
  prop       = allocate<par_prop,MEM_ALIGNMENT>(Parameter::SYS_SIZE);
  pv_bef     = allocate<double3,MEM_ALIGNMENT >(Parameter::SYS_SIZE);
  force	     = allocate<double3,MEM_ALIGNMENT >(Parameter::SYS_SIZE);
  force_bef  = allocate<double3,MEM_ALIGNMENT >(Parameter::SYS_SIZE);
  delta_sumr = allocate<double3,MEM_ALIGNMENT >(Parameter::SYS_SIZE);
#else
  pr         = new double3  [Parameter::SYS_SIZE];
  pv         = new double3  [Parameter::SYS_SIZE];
  prop       = new par_prop [Parameter::SYS_SIZE];
  pv_bef     = new double3  [Parameter::SYS_SIZE];
  force      = new double3  [Parameter::SYS_SIZE];
  force_bef  = new double3  [Parameter::SYS_SIZE];
  delta_sumr = new double3  [Parameter::SYS_SIZE];

#ifdef ADD_POSRES
  pr_base    = new double3  [Parameter::SYS_SIZE];
  rest_on    = new bool     [Parameter::SYS_SIZE];
#endif

#endif

  cheminfo.prtcl_chem	 = new bool    [Parameter::SYS_SIZE];
  cheminfo.lipid_chem	 = new bool    [param.ampN         ];
  cheminfo.lipid_idx	 = new int       [Parameter::SYS_SIZE];
  cheminfo.part_idx	 = new int       [Parameter::SYS_SIZE];
  cheminfo.lipid_unit	 = new int       [Parameter::SYS_SIZE];
  cheminfo.lip_elem_idx	 = new int       [param.ampN * Parameter::ALL_UNIT_N];
  cheminfo.head_elem_idx = new int       [Parameter::REAC_PART * param.headN];
  cheminfo.tail_elem_idx = new int       [Parameter::TAIL_PART * param.tailN];
  cheminfo.near_water	 = new float     [Parameter::SYS_SIZE];
  cheminfo.near_info	 = new float_int [Parameter::SYS_SIZE];

  celllist.near_prtcl_idx = new std::vector<int> [param.all_grid];
  for(int i = 0; i < param.all_grid; i++) celllist.near_prtcl_idx[i].resize(Parameter::BUF_SIZE);
  celllist.n_near_prtcl   = new int [param.all_grid];
  celllist.prtcl_idx      = new int [Parameter::SYS_SIZE];
  celllist.prtcl_in_cell  = new int [Parameter::SYS_SIZE];
}

void dpdsystem::DeallocMem() {
  delete p_observer;

#ifdef PROFILE_MODE
  delete p_profile;
#endif

  delete [] celllist.buck_elem;
  delete [] celllist.buck_addrs;
  delete [] celllist.next_dest;

#if 0
  deallocate<double3,MEM_ALIGNMENT>(pr);
  deallocate<double3,MEM_ALIGNMENT>(pv);
  deallocate<par_prop,MEM_ALIGNMENT>(prop);
  deallocate<double3,MEM_ALIGNMENT>(pv_bef);
  deallocate<double3,MEM_ALIGNMENT>(force);
  deallocate<double3,MEM_ALIGNMENT>(force_bef);
  deallocate<double3,MEM_ALIGNMENT>(delta_sumr);
#else  
  delete [] pr;
  delete [] pv;
  delete [] prop;
  delete [] pv_bef;
  delete [] force;
  delete [] force_bef;
  delete [] delta_sumr;

#ifdef ADD_POSRES
  delete [] pr_base;
  delete [] rest_on;
#endif

#endif

  delete [] cheminfo.prtcl_chem;
  delete [] cheminfo.lipid_chem;
  delete [] cheminfo.lipid_idx;
  delete [] cheminfo.part_idx;
  delete [] cheminfo.lipid_unit;
  delete [] cheminfo.lip_elem_idx;
  delete [] cheminfo.tail_elem_idx;
  delete [] cheminfo.head_elem_idx;
  delete [] cheminfo.near_water;
  delete [] cheminfo.near_info;

  delete [] celllist.near_prtcl_idx;
  delete [] celllist.n_near_prtcl;
  delete [] celllist.prtcl_idx;
  delete [] celllist.prtcl_in_cell;
}

void dpdsystem::RemoveCMDrift() {
  double3 cm_vel = std::accumulate(pv, pv + Parameter::SYS_SIZE, double3(0.0));
  cm_vel /= Parameter::SYS_SIZE;
  for (int i = 0; i < Parameter::SYS_SIZE; i++) pv[i] -= cm_vel;
}

#ifdef ADD_POSRES
void dpdsystem::SetPosresPrtcl() {
  for(int i = 0; i < Parameter::SYS_SIZE; i++) {
    if (!cheminfo.prtcl_chem[i] && cheminfo.lipid_unit[i] == (Parameter::ALL_UNIT_N - 2)) {
      rest_on[i] = true;
      pr_base[i] = pr[i];
    }
  }
}
#endif

void dpdsystem::CheckInitialized() const {
  const double3 zerov(0.0);
  for (int i = 0; i < Parameter::SYS_SIZE; i++) {
    CHECK_EQUATION(pr[i] <= p_param->L, pr[i]);
    CHECK_EQUATION(pr[i] >= zerov, pr[i]);
    CHECK_EQUATION(pv[i].isfinite3(), pv[i]);
    CHECK_EQUATION(prop[i] >= 0, prop[i]);
    CHECK_EQUATION(prop[i] < 3, prop[i]);
  }

  for (int i = 0; i < Parameter::SYS_SIZE; i++) {
    CHECK_EQUATION(cheminfo.lipid_idx[i] < p_param->ampN, cheminfo.lipid_idx[i]);
    CHECK_EQUATION((cheminfo.part_idx[i] < p_param->tailN) || (cheminfo.part_idx[i] < p_param->headN), cheminfo.part_idx[i]);
    CHECK_EQUATION(cheminfo.lipid_unit[i] < p_param->ALL_UNIT_N, cheminfo.lipid_unit[i]);
  }
  for (int i = 0; i < p_param->headN * Parameter::REAC_PART; i++)
    CHECK_EQUATION(cheminfo.head_elem_idx[i] < Parameter::SYS_SIZE, cheminfo.head_elem_idx[i]);
  for (int i = 0; i < p_param->tailN * Parameter::TAIL_PART; i++)
    CHECK_EQUATION(cheminfo.tail_elem_idx[i] < Parameter::SYS_SIZE, cheminfo.tail_elem_idx[i]);
  for (int i = 0; i < p_param->ampN * Parameter::ALL_UNIT_N; i++)
    CHECK_EQUATION(cheminfo.lip_elem_idx[i] < Parameter::SYS_SIZE, cheminfo.lip_elem_idx[i]);
  for (int i = 0; i < p_param->headN * Parameter::REAC_PART; i++)
    CHECK_EQUATION(prop[cheminfo.head_elem_idx[i]] != Water, prop[cheminfo.head_elem_idx[i]]);
  for (int i = 0; i < p_param->tailN * Parameter::TAIL_PART; i++) 
    CHECK_EQUATION(prop[cheminfo.tail_elem_idx[i]] != Water, prop[cheminfo.tail_elem_idx[i]]);
  for (int i = 0; i < p_param->ampN * Parameter::ALL_UNIT_N; i++) 
    CHECK_EQUATION(prop[cheminfo.lip_elem_idx[i]] != Water, prop[cheminfo.lip_elem_idx[i]]);
}

void dpdsystem::FirstStep(B_sorter& bsorter, ChemManager& chemmanage, F_calculator& f_calc, RNG& rng) {
  bsorter.MkPrtclIdx(pr, *p_param, celllist);
  bsorter.BucketSort(*this, cheminfo, *p_param, celllist);
  chemmanage.RegistLipidIdx(*p_param, cheminfo);
  bsorter.SetPrtclCellIdxNoSTL(celllist, *p_param);
  f_calc.AddConservForce(pr, prop, force, celllist, *p_param, cheminfo);
}

void dpdsystem::Execute(const int all_time, const int time_step_mic, const int time_step_mac, const int time_step_vt, const int chem_beg) {
  ChemManager chemmanage;
  F_calculator f_calc(*p_param);
  B_sorter bsorter(*p_param);
  RNG rng(1234, omp_get_max_threads());
  //RNG rng( (size_t)time(NULL), omp_get_max_threads() );

  p_initializer->GenParticles(*this, chemmanage, cheminfo, *p_param, rng);
#ifdef ADD_POSRES
  SetPosresPrtcl();
#endif

  CheckInitialized();
  
  FirstStep(bsorter, chemmanage, f_calc, rng);
  for (int time = 0; time < all_time; time++) {
    //std::cout << "time" <<  time << std::endl;
    if (time == Parameter::EQUIL_TIME) AdjustMomentum();

    VelocVerlet1();

    bsorter.MkPrtclIdx(pr, *p_param, celllist);
    if (time % B_sorter::SORT_FREQ == 0) {
      bsorter.BucketSort(*this,cheminfo,*p_param,celllist); 
      chemmanage.RegistLipidIdx(*p_param,cheminfo);
    }
    bsorter.SetPrtclCellIdxNoSTL(celllist,*p_param);


    f_calc.AddConservForce(pr, prop, force, celllist, *p_param, cheminfo);

#ifdef ADD_BIND_POTENT
    //if(time > 10000) f_calc.AddBindForceCyl(pr,prop,force,cheminfo,*p_param);
    //if(time > 0) f_calc.AddBindForceSph(pr,prop,force,cheminfo,*p_param);
    //if(time > 0) f_calc.AddBindForceZaxis(pr, prop, force, cheminfo, *p_param);
#endif

#ifdef ADD_POSRES
    if (time < 10000) f_calc.AddPosRes(pr, pr_base, force, rest_on);
#endif

#ifdef DEBUG
    bsorter.CheckSorted(*this,*p_param,celllist);
    bsorter.CheckNearList(*this,*p_param);
#endif

#ifdef CHEM_MODE
    if(time > chem_beg) chemmanage.ChemEvent(*p_param,cheminfo,rng);
    chemmanage.ClearNearInfo(*p_param,cheminfo);
#endif

    VelocVerlet2();

    if (time % Parameter::COL_FREQ == 0)
      f_calc.AddDisspRandom(pr, pv, prop, celllist, *p_param, rng);

    if (time % time_step_mic == 0) {
      p_observer->DumpTranject(*this, cheminfo, *p_param);
    }
    
    if (time % time_step_mac == 0) {
      p_observer->DumpMacroVal(*this, *p_param);
      p_observer->DumpEmbTailInfo(*this, cheminfo, *p_param, time);

#ifdef Y_REFLECT_BOUND
      p_observer->DumpLocalVal(*this, *p_param);
#endif
      
#ifdef CALC_HEIGHT
      p_observer->DumpMembHeight(*this, *p_param, time);
#endif

    }
    if (time % time_step_vt == 0) {
      p_observer->DumpPressure(*this, *p_param, Vir);
      p_observer->DumpConfigTempera(configT);
      
      UpdateObserverConfigT(f_calc.DumpConfigT(force) );
    }
  } //end of main loop
  
  //dump final configuration for restart simulation.
  p_observer->DumpFinalConfig(*this, cheminfo, *p_param);
}

