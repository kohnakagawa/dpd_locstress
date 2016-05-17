#pragma once

#include <iostream>
#include <iomanip>
#include <limits>
#include "mvector3.hpp"

struct PtclBuffer {
  double3 pos, vel;
  int prop = -1;
  bool p_chem = false, l_chem = false;
  int l_unit = -1, l_idx = -1, p_idx = -1;
  int hash = -1;

  static constexpr char atom_type[21] = {
    'O', 'N', 'C', 'S', 'P', 'Z', 'X', 'O', 'N', 'C', 'S', 'P', 'Z', 'X', 'O', 'N', 'C', 'S', 'P', 'Z', 'X'  
  };
  
  PtclBuffer() {
    pos = vel = double3(std::numeric_limits<double>::signaling_NaN());
  }
  
  void WriteOstream(std::ostream& ost) const {
    ost << std::setprecision(15);
    ost << atom_type[prop] << " " <<  pos << " " << vel
	<< " " << prop << " " << p_chem << " " << l_chem
	<< " " << l_unit << " " << l_idx << " " << p_idx << std::endl;
  }

  void ReadFromIstream(std::istream& ist) {
    char buf = '\0';
    ist >> buf >> pos.x >> pos.y >> pos.z >> vel.x >> vel.y >> vel.z >> prop >> p_chem
	>> l_chem >> l_unit >> l_idx >> p_idx;
  }
};
