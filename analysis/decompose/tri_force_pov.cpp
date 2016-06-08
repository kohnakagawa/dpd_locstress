#include "pov_renderer.hpp"
#include "ls_tools.hpp"
#include <random>
#include <iostream>

void gen_triangle(std::array<double3, 3>& pos,
		  std::array<double3, 3>& force,
		  std::array<double3, 3>& dF) {
  const double b_len = 4.0;
  
  pos[0].x = pos[0].y = pos[0].z = 0.0;
  pos[1].x = b_len; pos[1].y = pos[0].z = 0.0;
  pos[2].x = 0.5 * b_len; pos[2].y = 0.5 * b_len * std::sqrt(3.0);
  pos[2].z = 0.0;

  std::mt19937 mt(120);
  std::uniform_real_distribution<double> uni_dist_plus(0.0, 3.0);
  
  const double3 dr[] = {
    pos[1] - pos[0], // 0 -> 1
    pos[2] - pos[1], // 1 -> 2
    pos[0] - pos[2]  // 2 -> 0
  };

  std::array<double, 3> cf;

#if 0
  std::cout << "random coef\n";
  cf[0] = uni_dist_plus(mt);
  cf[1] = uni_dist_plus(mt);
  cf[2] = uni_dist_plus(mt);
#else
  cf.fill(1.5);
#endif

  std::cout << "force coef\n";
  std::cout << cf[0] << ", " << cf[1] << ", " << cf[2] << std::endl;
  
  dF[0] = cf[0] * dr[0];
  dF[1] = cf[1] * dr[1];
  dF[2] = cf[2] * dr[2];

  force[0] = dF[2] - dF[0];
  force[1] = dF[0] - dF[1];
  force[2] = dF[1] - dF[2];
}

int main(int argc, char* argv[]) {
  const std::string exe_modes[] {"sum_force", "pair_force", "center_force", "decompo_force"};
  const std::string file_names[] {"sum_force.pov", "pair_force.pov", "center_force.pov", "decompo_force.pov"};

  enum : int {
    SUM_FORCE = 0,
    PAIR_FORCE,
    CENTER_FORCE,
    DECOMP_FORCE,
      
    NUM_MODE
  };

  if (argc != 3) {
    std::cerr << "argv[1] == target directory name.\n";
    std::cerr << "argv[2] == execution mode.\n";
    for (const auto m : exe_modes)
      std::cout << m << " ";
    std::cout << "\n";
    std::exit(1);
  }

  const std::string mode = argv[2];

  std::array<double3, 3> pos, force, dF, dr;
  gen_triangle(pos, force, dF);
  
  const auto cm_pos = std::accumulate(pos.cbegin(),
				      pos.cend(),
				      double3 {0.0, 0.0, 0.0}) / pos.size();
  const std::string cur_dir = argv[1];
  int mode_id = -1;
  for (int i = 0; i < NUM_MODE; i++) if (mode == exe_modes[i]) mode_id = i;
  if (mode_id == -1) {
    std::cerr << "Please select correct execution mode.\n";
    for (const auto m : exe_modes) std::cout << m << " ";
    std::cout << "\n";
    std::exit(1);
  }
  const std::string fname = cur_dir + "/" + file_names[mode_id];
  
  PovRenderer<double3> renderer(fname.c_str());
  renderer.SetCamera(cm_pos + double3(0.0, 0.0, -24.0), cm_pos + double3(0.0, 0.0, 0.0));
  renderer.SetLight(cm_pos + double3(0.0, 0.0, -10.0));
  renderer.SetBackGround();
  renderer.AppendDefaultIncludeFiles();
  
  // atoms
  const double sphere_rad = 0.5;
  renderer.AppendObject(ObjectFactory<double3>::CreateSphere(pos[0], sphere_rad, PovColor::red<double3>()));
  renderer.AppendObject(ObjectFactory<double3>::CreateSphere(pos[1], sphere_rad, PovColor::red<double3>()));
  renderer.AppendObject(ObjectFactory<double3>::CreateSphere(pos[2], sphere_rad, PovColor::red<double3>()));
  
  // atom names
  const double offset = 1.0;
  const double scale  = 1.0;
  renderer.AppendObject(ObjectFactory<double3>::CreateText("1", "crystal.ttf", 0.10,
							   PovColor::black<double3>(), double3 {-0.6, 0.6, 0.0} * offset + pos[0], {scale, scale, scale}, {0.0, 0.0, 0.0}));
  renderer.AppendObject(ObjectFactory<double3>::CreateText("2", "crystal.ttf", 0.10,
							   PovColor::black<double3>(), double3 {0.6, 0.6, 0.0} * offset + pos[1], {scale, scale, scale}, {0.0, 0.0, 0.0}));
  renderer.AppendObject(ObjectFactory<double3>::CreateText("3", "crystal.ttf", 0.10,
							   PovColor::black<double3>(), double3 {0.8, -0.8, 0.0} * offset + pos[2], {scale, scale, scale}, {0.0, 0.0, 0.0}));
  
  dr[0] = pos[1] - pos[0];
  dr[1] = pos[2] - pos[1];
  dr[2] = pos[0] - pos[2];
  
  // offset for phantom bonds
  const double3 offset_vec(0.0, 0.0, 0.2);

  // add force vectors
  const double f_scale = 0.20;
  const double cyl_rad = 0.085;

  if (mode == exe_modes[SUM_FORCE]) {
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[0] + offset_vec, dr[0], {1.0, 1.0, 1.0}, "M_Glass", 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[1] + offset_vec, dr[1], {1.0, 1.0, 1.0}, "M_Glass", 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[2] + offset_vec, dr[2], {1.0, 1.0, 1.0}, "M_Glass", 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[0],
    							     force[0] * f_scale,
    							     cyl_rad, PovColor::black<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[1],
    							     force[1] * f_scale,
    							     cyl_rad, PovColor::black<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[2],
    							     force[2] * f_scale,
    							     cyl_rad, PovColor::black<double3>()));
  } else if (mode == exe_modes[PAIR_FORCE]) {
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[0] + offset_vec, dr[0], {1.0, 1.0, 1.0}, "M_Glass", 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[1] + offset_vec, dr[1], {1.0, 1.0, 1.0}, "M_Glass", 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[2] + offset_vec, dr[2], {1.0, 1.0, 1.0}, "M_Glass", 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[0],
    							     -dF[0] * f_scale,
    							     cyl_rad, PovColor::black<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[0],
    							     dF[2] * f_scale,
    							     cyl_rad, PovColor::blue<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[1],
    							     dF[0] * f_scale,
    							     cyl_rad, PovColor::black<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[1],
    							     -dF[1] * f_scale,
    							     cyl_rad, PovColor::green<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[2],
    							     dF[1] * f_scale,
    							     cyl_rad, PovColor::green<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[2],
    							     -dF[2] * f_scale,
    							     cyl_rad, PovColor::blue<double3>()));
  } else if (mode == exe_modes[CENTER_FORCE]) {
    renderer.AppendObject(ObjectFactory<double3>::CreateSphere(cm_pos, sphere_rad, PovColor::green<double3>()));
    
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(cm_pos + offset_vec, pos[0] - cm_pos, {1.0, 1.0, 1.0}, "M_Glass", 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(cm_pos + offset_vec, pos[1] - cm_pos, {1.0, 1.0, 1.0}, "M_Glass", 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(cm_pos + offset_vec, pos[2] - cm_pos, {1.0, 1.0, 1.0}, "M_Glass", 0.1));
    
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[0],
							       force[0] * f_scale,
							       cyl_rad, PovColor::black<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(cm_pos,
							       -force[0] * f_scale,
							       cyl_rad, PovColor::black<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[1],
							       force[1] * f_scale,
							       cyl_rad, PovColor::blue<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(cm_pos,
							       -force[1] * f_scale,
							       cyl_rad, PovColor::blue<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[2],
							       force[2] * f_scale,
							       cyl_rad, PovColor::yellow<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(cm_pos,
							       -force[2] * f_scale,
							       cyl_rad, PovColor::yellow<double3>()));
  } else if (mode == exe_modes[DECOMP_FORCE]) {
    double f_decomp_err = 1.0;
    std::array<double3, 5> dF_d;
    const double3 base_pos = ((pos[1] + pos[2]) * 0.5) * 0.5 + ((pos[0] + pos[1]) * 0.5) * 0.5;
    decompose_3n(pos[0], pos[1], pos[2], base_pos, force[0], force[1], force[2], f_decomp_err, &dF_d[0]);

    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[0] + offset_vec, base_pos - pos[0], {1.0, 1.0, 1.0}, "M_Glass", 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[1] + offset_vec, base_pos - pos[1], {1.0, 1.0, 1.0}, "M_Glass", 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[2] + offset_vec, base_pos - pos[2], {1.0, 1.0, 1.0}, "M_Glass", 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[2] + offset_vec, pos[1] - pos[2], {1.0, 1.0, 1.0}, "M_Glass", 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[1] + offset_vec, pos[0] - pos[1], {1.0, 1.0, 1.0}, "M_Glass", 0.1));

    renderer.AppendObject(ObjectFactory<double3>::CreateSphere(base_pos, sphere_rad, PovColor::green<double3>()));
  
    // atom 1
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[0],
    							       dF_d[0] * f_scale,
    							       cyl_rad, PovColor::blue<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[0],
    							       dF_d[1] * f_scale,
    							       cyl_rad, PovColor::black<double3>()));
    // atom 2
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[1],
    							       -dF_d[0] * f_scale,
    							       cyl_rad, PovColor::blue<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[1],
    							       dF_d[2] * f_scale,
    							       cyl_rad, PovColor::yellow<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[1],
    							       dF_d[3] * f_scale,
    							       cyl_rad, PovColor::cyan<double3>()));
    // atom 3
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[2],
    							       -dF_d[2] * f_scale,
    							       cyl_rad, PovColor::yellow<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[2],
    							       dF_d[4] * f_scale,
    							       cyl_rad, PovColor::magenta<double3>()));
    // atom 4
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(base_pos,
    							       -dF_d[1] * f_scale,
    							       cyl_rad, PovColor::black<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(base_pos,
    							       -dF_d[3] * f_scale,
    							       cyl_rad, PovColor::cyan<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(base_pos,
    							       -dF_d[4] * f_scale,
    							       cyl_rad, PovColor::magenta<double3>()));
  }

  renderer.WriteRenderScript();
}
