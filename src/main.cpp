#include "dpdsystem.hpp"

#ifdef DEBUG
#warning "DEBUG MODE IS SET!"
#endif

#ifdef PROFILE_MODE
#warning "PROFILE MODE IS SET!"
#endif

#ifndef CHEM_MODE
#warning "THERE IS NO CHEMICAL REACTION!"
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
  const int all_time = 100, time_step_mic = 100, time_step_mac = 100, time_step_vt = 100, chem_beg = 1000;
  dpdsys.Execute(all_time, time_step_mic, time_step_mac, time_step_vt, chem_beg);

  const std::string fname = std::string(argv[1]) + "/all_param.txt";
  std::ofstream fout(fname.c_str());
  param.DumpAllParam(fout);
  param.DumpGraphicDat(all_time, time_step_mic);
}
