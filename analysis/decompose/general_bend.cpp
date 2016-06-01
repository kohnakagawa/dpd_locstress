#include "ls_tools.hpp"
#include <random>
#include <ctime>
#include <iomanip>
#include <map>

void gen_init_config(std::array<double3, 3>& pos,
		     std::array<double3, 3>& force,
		     const std::array<double, 2> b_len,
		     const double theta) {
  pos[0].x = pos[0].y = pos[0].z = 0.0;
  pos[1] = pos[0];
  pos[1].x += b_len[0];
  pos[2] = pos[1] + double3(b_len[1] * std::cos(theta), b_len[1] * std::sin(theta), 0.0);
  
  std::mt19937 mt(static_cast<size_t>(time(nullptr)));
  // std::uniform_real_distribution<double> uni_dist(-10.0, 10.0);
  std::uniform_real_distribution<double> uni_dist_plus(0.0, 10.0);
  std::uniform_real_distribution<double> uni_dist_minus(-10.0, 0.0);

  std::cout << std::setprecision(15);

  const double3 dr[] = {
    pos[1] - pos[0], // 0 -> 1
    pos[2] - pos[1], // 1 -> 2
    pos[0] - pos[2]  // 2 -> 0
  };

  double3 dF[3];
  std::cout << "random coef\n";
  std::array<double, 3> cf;
  cf[0] = uni_dist_plus(mt);
  cf[1] = uni_dist_minus(mt);
  cf[2] = uni_dist_plus(mt);

  std::cout << "force coef\n";
  std::cout << cf[0] << ", " << cf[1] << ", " << cf[2] << std::endl;
  
  dF[0] = cf[0] * dr[0];
  dF[1] = cf[1] * dr[1];
  dF[2] = cf[2] * dr[2];
  
  // std::cout << "ana solution: " << calc_stress_sum_cent(cf, dr) << std::endl;
  const auto min_pos = get_stress_sum_min_pos(dr, pos[1]);
  
  const double3 F[] {
    -dF[0] + dF[2],
     dF[0] - dF[1],
     dF[1] - dF[2]};
  
  double3 dF_d[5];
  double f_decomp_err = 0.0;
  decompose_3n(pos[0], pos[1], pos[2], min_pos, F[0], F[1], F[2], f_decomp_err, dF_d);
  
  double3 dr_new[] {
    pos[0] - pos[1],
    pos[0] - min_pos,
    pos[1] - pos[2],
    pos[1] - min_pos,
    pos[2] - min_pos
  };

  std::cout << "analysis at minimum stress position.\n";
  std::cout << "|dF14||dr14|, |dF24||dr24|, |dF34||dr34|\n";  
  std::cout << dF_d[1].norm2() * dr_new[1].norm2() << ", " << dF_d[4].norm2() * dr_new[4].norm2() << ", " << dF_d[3].norm2() * dr_new[3].norm2() << std::endl;
  std::cout << "|dF12||dr12|, |dF23||dr23|\n";
  std::cout << dF_d[0].norm2() * dr_new[0].norm2() << " " << dF_d[2].norm2() * dr_new[2].norm2() << std::endl;
  std::cout << "end of analysis\n";
  
  force[0] = dF[2] - dF[0];
  force[1] = dF[0] - dF[1];
  force[2] = dF[1] - dF[2];

  std::cout << "CHECK of 3 body potential.\n";
  std::cout << "linear momentum " << force[0] + force[1] + force[2] << std::endl;
  std::cout << "angular momentum " << (force[0] ^ pos[0]) + (force[1] ^ pos[1]) + (force[2] ^ pos[2]) << std::endl;

  std::cout << "stress sum (central_force/dipole) ";
  std::cout << calc_stress_sum(dF, dr, 3) << " ";
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "argv[1] == target directory name.\n";
    std::exit(1);
  }

  std::string cur_dir = argv[1];
  
  LSParam lsparam;
  lsparam.LoadParam(argv[1]);
  
  std::array<double3, 3> pos, force;
  gen_init_config(pos, force, lsparam.bleng, lsparam.bend_angle);
  
  double3 dF[] = {
    double3(0.0, 0.0, 0.0), // 2 -> 1
    double3(0.0, 0.0, 0.0), // 4 -> 1
    double3(0.0, 0.0, 0.0), // 3 -> 2
    double3(0.0, 0.0, 0.0), // 4 -> 2
    double3(0.0, 0.0, 0.0)  // 4 -> 3
  };
  
  double3 dr_new[] = {
    double3(0.0, 0.0, 0.0), // 2 -> 1
    double3(0.0, 0.0, 0.0), // 4 -> 1
    double3(0.0, 0.0, 0.0), // 3 -> 2
    double3(0.0, 0.0, 0.0), // 4 -> 2
    double3(0.0, 0.0, 0.0)  // 4 -> 3
  };

  dr_new[0] = pos[0] - pos[1];
  dr_new[2] = pos[1] - pos[2];
  
  const auto relaxed_pos = get_relaxed_pos(pos[0], pos[1], force[0], force[1]);
  
  dr_new[1] = pos[0] - relaxed_pos;
  dr_new[3] = pos[1] - relaxed_pos;
  dr_new[4] = pos[2] - relaxed_pos;  

  double f_decomp_err = 0.0;
  decompose_3n(pos[0], pos[1], pos[2], relaxed_pos, force[0], force[1], force[2], f_decomp_err, dF);

  std::cout << calc_stress_sum(dF, dr_new, 5) << std::endl;

  const double grid_sp = 0.037;
  const double3 dw_pos(-0.9, -0.9, 0.0);
  const double3 up_pos = pos[2] + double3(0.9, 0.9, 0.0);
  std::array<int, 2> idx = {{0, 0}};

  std::string fname = cur_dir + "/Fr_stress_dist.txt";
  std::ofstream fout(fname.c_str());
  fout << std::setprecision(10);

  std::pair<double, double3> smin_w_pos;
  smin_w_pos.first = std::numeric_limits<double>::max();
  while (true) {
    const auto base_pos = dw_pos + double3(idx[0] * grid_sp, idx[1] * grid_sp, 0.0);
    idx[0]++;
    if (base_pos[0] > up_pos[0]) {
      idx[0] = 0;
      idx[1]++;
    }

    if (base_pos[1] > up_pos[1]) {
      break;
    }

    decompose_3n(pos[0], pos[1], pos[2], base_pos, force[0], force[1], force[2], f_decomp_err, dF);
    
    dr_new[1] = pos[0] - base_pos;
    dr_new[3] = pos[1] - base_pos;
    dr_new[4] = pos[2] - base_pos;

    const auto stress_sum = calc_stress_sum(dF, dr_new, 5);
    fout << base_pos[0] << " " << base_pos[1] << " " << stress_sum << std::endl;
    if (smin_w_pos.first > stress_sum) {
      smin_w_pos.first = stress_sum;
      smin_w_pos.second = base_pos;
    }

    if (base_pos[0] > up_pos[0]) fout << std::endl;
  }
  
  fout.close();
  fname = cur_dir + "/Fr_cent.txt";
  fout.open(fname.c_str());
  fout << std::setprecision(10);
  const auto dr_cent = relaxed_pos - pos[1];
  fout << pos[1] << " " << dr_cent << std::endl;
  
  const auto bi_vec = (pos[0] + pos[2]) * 0.5 - pos[1];
  fout << pos[1] << " " << bi_vec << std::endl;
  
  const auto min_vec = smin_w_pos.second - pos[1];
  fout << pos[1] << " " << min_vec << " #" << smin_w_pos.first << std::endl;

  fout.close();
  fname = cur_dir + "/stress_min_pos.txt";
  fout.open(fname.c_str(), std::ios::app);
  fout << smin_w_pos.second << " #";
}
