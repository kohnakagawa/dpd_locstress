#include "dpdsystem.hpp"

#ifdef DEBUG
#warning "DEBUG MODE IS SET!"
#endif

#ifdef PROFILE_MODE
#warning "PROFILE MODE IS SET!"
#endif

#ifndef CHEM_MODE
#warning "THERE IS NO CHEMICAL REACTION!"
#else
#error "THIS MODE IS DEPRECATED IN CURRENT VERSION!"
#endif

#ifdef Y_REFLECT_BOUND
#warning "BOUNDARY CONDITION IS NOT PERIODIC BOUNDARY ALONG Y AXIS."
#endif

#ifdef NO_THERMO_STAT
#warning "THERE IS NO DPD THERMOSTAT."
#endif

#ifdef VISC_KERN_SQRT
#warning "DISSIPATION KERNEL FUNCION IS SQRT."
#endif

#ifdef ADD_BIND_POTENT
#warning "BIND POTENTIAL IS ADDED."
#endif

#ifdef CALC_LOC_STRESS
#warning "WILL CALCULATE LOCAL STRESS."
#endif

namespace {
  void warning(const int argc) {
    if (argc != 2) {
      std::cerr << "Usage:\n";
      std::cerr << "argv[1] = target directory\n";
      std::exit(1);
    }
  }
}

int main(int argc, char *argv[]) {
  warning(argc);

  Parameter param(argv[1]);
  param.LoadParam();
  param.LoadCheck();

  dpdsystem dpdsys(param);
  dpdsys.Initialize();
  dpdsys.Execute(param.all_time(), param.time_step_mic(), param.time_step_mac(), param.time_step_vt(), param.chem_beg());

  const std::string fname = std::string(argv[1]) + "/all_param.txt";
  std::ofstream fout(fname.c_str());
  param.DumpAllParam(fout);
  param.DumpGraphicDat(param.all_time(), param.time_step_mic());
}
