#include "../../src/force_calculator.hpp"
#include <iostream>
#include <memory>

int Parameter::sys_size;

void gen_init_config(std::array<double3, 3>& pos,
		     const double b_len,
		     const double theta) {
  pos[0].x = pos[0].y = pos[0].z = 0.0;
  pos[1] = pos[0];
  pos[1].x += b_len;
  pos[2] = pos[1] + double3(b_len * std::cos(theta), b_len * std::sin(theta), 0.0);
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
    double3(0.0, 0.0, 0.0),
    double3(0.0, 0.0, 0.0)
  };
  
  F_calculator fcalc(param);
  
  tensor3d d_virial = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double   d_lap = 0.0;
  fcalc.StoreBendForce(&pos[0], dr, inv_dr, dist, d_virial, d_lap, F, param);
  std::cout << fcalc.DumpFdecompError() << std::endl;
  
  const auto loc_stress = fcalc.DumpCurLocStress();
  const int thnum = loc_stress.size();
  const int grid_num = loc_stress[0].size();
  std::unique_ptr<tensor3d[]> loc_stress_sum(new tensor3d [grid_num]);
  for (int i = 0; i < thnum; i++)
    for (int g = 0; g < grid_num; g++)
      for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) {
	  loc_stress_sum[g][j][k] += loc_stress[i][g][j][k];
	}
  
  tensor3d all_stress = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int g = 0; g < grid_num; g++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) 
	all_stress[j][k] += loc_stress_sum[g][j][k];

  std::cout << "CHECK: local stress sum, global stress\n";
  for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) {
      std::cout << all_stress[j][k] << ", " << d_virial[j][k] << std::endl;
    }
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "argv[1] == target directory name.\n";
    std::exit(1);
  }
  
  Parameter param(argv[1]);
  param.LoadParam();
  param.LoadCheck();
  
  std::array<double3, 3> pos;
  gen_init_config(pos, 1.0, M_PI / 6.0);
  
  calc_force(pos, param);
}
