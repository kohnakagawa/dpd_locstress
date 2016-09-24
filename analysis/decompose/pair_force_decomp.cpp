#include "pov_renderer.hpp"
#include "ls_tools.hpp"
#include <random>
#include <iostream>

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "argv[1] == target directory name.\n";
    std::cerr << "argv[2] == original or decomp.\n";
    std::exit(1);
  }

  const std::string cur_dir = argv[1];
  const std::string exe_mode = argv[2];

  if (exe_mode != "original" && exe_mode != "decomp") {
    std::cerr << "execution mode should be original or decomp.\n";
    std::exit(1);
  }

  // run input paramter
  const double b_len = 6.0, angle = 40.0 * M_PI / 180.0, cf = 1.5;

  std::array<double3, 4> pos;
  pos[0] = double3(0.0, 0.0, 0.0);
  pos[1] = double3(b_len, 0.0, 0.0);
  pos[2] = double3(0.5 * b_len, 0.5 * b_len * std::tan(angle), 0.0);
  pos[3] = double3(0.5 * b_len, -0.5 * b_len * std::tan(angle), 0.0);

  std::array<double3, 5> dF;
  dF[0] = cf * 0.5 * double3(1.0, -std::tan(angle), 0.0);
  dF[1] = cf * 0.5 * double3(1.0, std::tan(angle), 0.0);
  dF[2] = cf * 0.5 * double3(-1.0, -std::tan(angle), 0.0);
  dF[3] = cf * 0.5 * double3(-1.0, std::tan(angle), 0.0);
  dF[4] = cf * double3(0.0, std::tan(angle), 0.0);

  std::array<double3, 5> dr;
  dr[0] = pos[2] - pos[0];
  dr[1] = pos[3] - pos[0];
  dr[2] = pos[2] - pos[1];
  dr[3] = pos[3] - pos[1];
  dr[4] = pos[3] - pos[2];

  const std::string fname = cur_dir + "/pair_decomp.pov";
  PovRenderer<double3> renderer(fname.c_str());
  const auto cm_pos = std::accumulate(pos.cbegin(), pos.cend(), double3(0.0, 0.0, 0.0)) / pos.size();
  renderer.SetCamera(cm_pos + double3(0.0, 0.0, -24.0), cm_pos + double3(0.0, 0.0, 0.0));
  renderer.SetLight(cm_pos + double3(0.0, 0.0, -10.0));
  renderer.SetBackGround();
  renderer.AppendDefaultIncludeFiles();

  // atoms
  const double sphere_rad = 0.5;
  renderer.AppendObject(ObjectFactory<double3>::CreateSphere(pos[0],
                                                             sphere_rad,
                                                             PovColor::brown<double3>()));
  renderer.AppendObject(ObjectFactory<double3>::CreateSphere(pos[1],
                                                             sphere_rad,
                                                             PovColor::brown<double3>()));
  renderer.AppendObject(ObjectFactory<double3>::CreateSphere(pos[2],
                                                             sphere_rad,
                                                             PovColor::green<double3>()));
  renderer.AppendObject(ObjectFactory<double3>::CreateSphere(pos[3],
                                                             sphere_rad,
                                                             PovColor::green<double3>()));

  // bonds
  const double3 offset_vec(0.0, 0.0, 0.2);
  const double bond_rad = 0.1;
  if (exe_mode == "original") {
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[0] + offset_vec,
                                                                 pos[1] - pos[0],
                                                                 {1.0, 1.0, 1.0},
                                                                 "M_Glass",
                                                                 0.1));
  } else if (exe_mode == "decomp") {
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[0] + offset_vec,
                                                                 dr[0],
                                                                 {1.0, 1.0, 1.0},
                                                                 "M_Glass",
                                                                 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[0] + offset_vec,
                                                                 dr[1],
                                                                 {1.0, 1.0, 1.0},
                                                                 "M_Glass",
                                                                 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[1] + offset_vec,
                                                                 dr[2],
                                                                 {1.0, 1.0, 1.0},
                                                                 "M_Glass",
                                                                 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[1] + offset_vec,
                                                                 dr[3],
                                                                 {1.0, 1.0, 1.0},
                                                                 "M_Glass",
                                                                 0.1));
    renderer.AppendObject(ObjectFactory<double3>::CreateCylinder(pos[2] + offset_vec,
                                                                 dr[4],
                                                                 {1.0, 1.0, 1.0},
                                                                 "M_Glass",
                                                                 0.1));
  }

  // force vectors
  const double f_scale = 1.8;
  const double cyl_rad = 0.085;
  if (exe_mode == "original") {
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[0],
                                                               (dF[0] + dF[1]) * (-f_scale),
                                                               cyl_rad,
                                                               PovColor::black<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[1],
                                                               (dF[2] + dF[3]) * (-f_scale),
                                                               cyl_rad,
                                                               PovColor::black<double3>()));
  } else if (exe_mode == "decomp") {
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[3],
                                                               dF[0] * f_scale,
                                                               cyl_rad,
                                                               PovColor::blue<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[0],
                                                               -dF[0] * f_scale,
                                                               cyl_rad,
                                                               PovColor::blue<double3>()));

    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[3],
                                                               dF[2] * f_scale,
                                                               cyl_rad,
                                                               PovColor::red<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[1],
                                                               -dF[2] * f_scale,
                                                               cyl_rad,
                                                               PovColor::red<double3>()));

    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[2],
                                                               dF[1] * f_scale,
                                                               cyl_rad,
                                                               PovColor::black<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[0],
                                                               -dF[1] * f_scale,
                                                               cyl_rad,
                                                               PovColor::black<double3>()));

    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[2],
                                                               dF[3] * f_scale,
                                                               cyl_rad,
                                                               PovColor::cyan<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[1],
                                                               -dF[3] * f_scale,
                                                               cyl_rad,
                                                               PovColor::cyan<double3>()));

    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[3],
                                                               dF[4] * f_scale,
                                                               cyl_rad,
                                                               PovColor::orange_yellow<double3>()));
    renderer.AppendObject(ObjectFactory<double3>::CreateVector(pos[2],
                                                               -dF[4] * f_scale,
                                                               cyl_rad,
                                                               PovColor::orange_yellow<double3>()));
  }

  renderer.WriteRenderScript();
}
