#include "../../src/force_calculator.hpp"
#include <iostream>

void gen_init_config(std::array<double3, 3>& pos,
		     std::array<double3, 3>& force,
		     const double b_len,
		     const double theta) {
  pos[0].x = pos[0].y = pos[0].z = 0.0;
  pos[1] = pos[0];
  pos[1].y += b_len;
  pos[2] = pos[1] + double3(b_len * std::cos(theta), b_len * std::sin(theta), 0.0);
  force.fill(double3(0.0, 0.0, 0.0));
}

void calc_force(std::array<double3, 3>& pos,
		const Parameter& param) {
  const double3 dr[] = {
    pos[1] - pos[0],
    pos[2] - pos[1]
  };
  
  const double dist[] = {
    dr[0].dist2(),
    dr[1].dist2()
  };
  
  const double inv_dr[] = {
    1.0 / std::sqrt(dist[0]),
    1.0 / std::sqrt(dist[1])
  };
  
  double3 F[] = {
    double3(0.0, 0.0, 0.0),
    double3(0.0, 0.0, 0.0)
  };
  
  F_calculator fcalc(param);
  
  tensor3d d_virial = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double   d_lap = 0.0;
  fcalc.StoreBendForce(dr, inv_dr, dist, d_virial, d_lap, F);
  
  
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "argv[1] == target directory name.\n";
    std::exit(1);
  }
  
  Parameter param(argv[1]);
  param.LoadParam();
  param.LoadCheck();
  
  std::array<double3, 3> pos, force;
  gen_init_config(pos, force, 1.0, M_PI / 6.0);
  
  
}

