#include "ls_tools.hpp"
#include <iostream>
#include <memory>
#include <cstdlib>
#include <sstream>
#include <numeric>
#include <limits>
#include <iomanip>
#include <random>

int Parameter::sys_size;

const char modes[][256] = {
  "mspos",
  "bisect",
  "allof"
};

void gen_init_config(std::array<double3, 3>& pos,
		     const std::array<double, 2> b_len,
		     const double theta) {
  pos[0].x = pos[0].y = pos[0].z = 0.0;
  pos[1] = pos[0];
  pos[1].x += b_len[0];
  pos[2] = pos[1] + double3(b_len[1] * std::cos(theta), b_len[1] * std::sin(theta), 0.0);

  std::cout << "particle positions are\n";
  for (int i = 0; i < 3; i++)
    std::cout << "(" << pos[i] << "), ";
  std::cout << "\n\n";
}

void calc_cosine_bend(const double3* r,
		      const double3* dr,
		      const double*  inv_dr,
		      const double3& base_pos,
		      tensor3d& d_virial,
		      double3* dF,
		      double3* sumF,
		      double& f_decomp_err,
		      const Parameter& param) {
  const double	inv_dr_prod = inv_dr[0] * inv_dr[1];
  const double	inv_dist[2] = {inv_dr[0] * inv_dr[0],
			       inv_dr[1] * inv_dr[1]};
  const double	in_prod	    = dr[0] * dr[1];
  const double	cf_b	    = param.GetIntractions().cf_bend * inv_dr_prod;
  const double        cf_crs[2]   = {in_prod * inv_dist[0],
				     in_prod * inv_dist[1]};

  const double3 Ftb0(cf_b * (dr[1].x - cf_crs[0] * dr[0].x),
		     cf_b * (dr[1].y - cf_crs[0] * dr[0].y),
		     cf_b * (dr[1].z - cf_crs[0] * dr[0].z));
  const double3 Ftb1(cf_b * (dr[0].x - cf_crs[1] * dr[1].x),
		     cf_b * (dr[0].y - cf_crs[1] * dr[1].y),
		     cf_b * (dr[0].z - cf_crs[1] * dr[1].z));
    
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
      d_virial[i][j] += dr[0][i] * Ftb0[j] + dr[1][i] * Ftb1[j];	
    }
    
  const double3 Ftb_sum = Ftb0 - Ftb1;
  
  decompose_3n(r[0], r[1], r[2], base_pos, -Ftb0, Ftb_sum, Ftb1, f_decomp_err, dF);

  sumF[0] = -Ftb0;
  sumF[1] = Ftb_sum;
  sumF[2] = Ftb1;
}

void analyze_at_minimum_stress_pos(const double3* pos,
				   const double3* dr,
				   const double* inv_dr,
				   const double3& ms_pos,
				   double3* dr_new,
				   tensor3d& d_virial,
				   double3* dF,
				   double3* sumF,
				   double& f_decomp_err,
				   const Parameter& param) {
  calc_cosine_bend(pos, dr, inv_dr, ms_pos, d_virial, dF, sumF, f_decomp_err, param);

  dr_new[1] = pos[0] - ms_pos;
  dr_new[3] = pos[1] - ms_pos;
  dr_new[4] = pos[2] - ms_pos;
  
  std::cout << "analysis at minimum stress position.\n";
  std::cout << "|dF14||dr14|, |dF24||dr24|, |dF34||dr34|\n";
  std::cout << dF[1].norm2() * dr_new[1].norm2() << ", "
	    << dF[3].norm2() * dr_new[3].norm2() << ", "
	    << dF[4].norm2() * dr_new[4].norm2() << "\n\n";
  
  const auto cf = dr_new[1].norm2() / dr_new[4].norm2();
  std::cout << "length cf is " << cf << std::endl;

  const auto dF24_norm_1 = cf * dF[0].norm2() * (dF[0] * dr_new[0]) / (dF[0].norm2() * dr_new[0].norm2());
  const auto dF24_norm_2 = dF[2].norm2() / cf * (dF[2] * dr_new[2]) / (dF[2].norm2() * dr_new[2].norm2());

  std::cout << "decompose dF24 force.\n";
  
  std::cout << "dF24_norm_1 " << dF24_norm_1 << std::endl;
  std::cout << "dF24_norm_2 " << dF24_norm_2 << std::endl;
  
  std::cout << "original dF24 norm, reconstructed norm.\n";
  std::cout << dF[3].norm2() << ", " << std::abs(dF24_norm_1 + dF24_norm_2) << std::endl;

  // std::cout << dF[0].norm2() * dr_new[0].norm2() << " "<< dF24_norm_1 * dr_new[3].norm2() << std::endl;
  // std::cout << dF24_norm_2 * dr_new[3].norm2() << " " << dF[2].norm2() * dr_new[2].norm2() << std::endl;
  // std::cout << dF[0].norm2() * dr_new[0].norm2() << " " << dF[2].norm2() * dr_new[2].norm2() << std::endl;

  check_congruence_of_two_tri(pos, dr[0] * dr[0], dr[1] * dr[1], ms_pos);

  std::cout << "end of analysis.\n";
}


void show_result(const double3* sumF,
		 const double3* dF,
		 const double f_decomp_err,
		 const double3* dr,
		 const double lambda,
		 double& lambda_at_lowest_stress,
		 std::ofstream& fout) {
  static double lowest_stress = std::numeric_limits<double>::max();
  
  fout << lambda << " ";
  
  for (int i = 0; i < 5; i++) {
    fout << dF[i] * (dr[i] / dr[i].norm2()) << " ";
  }
  
  const auto sigma_sum = calc_stress_sum(dF, dr, 5);
  fout << sigma_sum << " ";

  if (lambda > 0.0) {
    if (lowest_stress > sigma_sum) {
      lowest_stress = sigma_sum;
      lambda_at_lowest_stress = lambda;
    }
  }
  
  const auto dx_new1 = -sumF[0] / sumF[0].norm2();
  const auto dy_new1 = dF[0] / dF[0].norm2();
  
  const auto dx_new2 = sumF[2] / sumF[2].norm2();
  const auto dy_new2 = -dF[2] / dF[2].norm2();

  const double sig_xx_12 = (dF[0] * dx_new1) * (dr[0] * dx_new1) + (dF[1] * dx_new1) * (dr[1] * dx_new1) + 0.5 * (dF[3] * dx_new1) * (dr[3] * dx_new1);
  const double sig_yy_12 = (dF[0] * dy_new1) * (dr[0] * dy_new1) + (dF[1] * dy_new1) * (dr[1] * dy_new1) + 0.5 * (dF[3] * dy_new1) * (dr[3] * dy_new1);
  
  const double sig_xx_23 = (dF[2] * dx_new2) * (dr[2] * dx_new2) + 0.5 * (dF[3] * dx_new2) * (dr[3] * dx_new2) + (dF[4] * dx_new2) * (dr[4] * dx_new2);
  const double sig_yy_23 = (dF[2] * dy_new2) * (dr[2] * dy_new2) + 0.5 * (dF[3] * dy_new2) * (dr[3] * dy_new2) + (dF[4] * dy_new2) * (dr[4] * dy_new2);

  fout << sig_xx_12 << " " << sig_yy_12 << " " << sig_xx_23 << " " << sig_yy_23 << " ";
  fout << dr[4] << " ";

  fout << f_decomp_err;
}

void show_result(const double3* dF,
		 const double3* dr,
		 const double3& base_pos,
		 std::ofstream& fout) {
  const auto sigma_sum = calc_stress_sum(dF, dr, 5);
  fout << base_pos[0] << " " << base_pos[1] << " " << sigma_sum << "\n";
}

void analyze_along_bisect(const double3* pos,
			  const double3* dr,
			  const double* inv_dr,
			  const double3& ms_pos,
			  double3* dr_new,
			  tensor3d& d_virial,
			  double3* dF,
			  double3* sumF,
			  double& f_decomp_err,
			  const Parameter& param,
			  std::ofstream& fout) {
  fout << std::setprecision(10);
  fout << "# lambda \t|dF21| \t|dF41| \t|dF32| \t|dF42| \t|dF43| \t|dF||dr| \tsigma_xx_21 \tsigma_yy_21 \tsigma_xx_23 \tsigma_yy_23 \tf_decomp_err \t dr43 \t |dF41||dr41|/(|dF43|/|dr43|)\n";

  double lambda_at_lowest_stress = 0.0;
  double3 base_vec = ms_pos - pos[1];

  for (double lam = -4.0; lam < 4.0; lam += 0.0015) {
    const double3 base_pos = pos[1] + base_vec * lam;
    calc_cosine_bend(&pos[0], dr, inv_dr, base_pos, d_virial, dF, sumF, f_decomp_err, param);
    
    dr_new[1] = pos[0] - base_pos;
    dr_new[3] = pos[1] - base_pos;
    dr_new[4] = pos[2] - base_pos;
    
    show_result(sumF, dF, f_decomp_err, dr_new, lam, lambda_at_lowest_stress, fout);
    fout << dF[1].norm2() * dr_new[1].norm2() / (dF[4].norm2() * dr_new[4].norm2()) << std::endl;
  }

  std::cout << std::setprecision(10);
  std::cout << "lambda at lowest stress = " << lambda_at_lowest_stress << std::endl;
}

void analyze_all_of_plane(const double3* pos,
			  const double3* dr,
			  const double* inv_dr,
			  double3* dr_new,
			  tensor3d& d_virial,
			  double3* dF,
			  double3* sumF,
			  double& f_decomp_err,
			  const Parameter& param,
			  std::ofstream& fout) {
  fout << std::setprecision(10);
  
  const double grid_sp = 0.037;
  
  const double3 dw_pos(-0.9, -0.9, 0.0);
  const double3 up_pos = pos[2] + double3(3.0, 3.0, 0.0);
  
  std::array<int, 2> idx = {{0, 0}};
  while (true) {
    const double3 base_pos = dw_pos + double3(idx[0] * grid_sp, idx[1] * grid_sp, 0.0);
    idx[0]++;
    if (base_pos[0] > up_pos[0]) {
      idx[0] = 0;
      idx[1]++;
    }

    if (base_pos[1] > up_pos[1]) {
      break;
    }
    
    calc_cosine_bend(&pos[0], dr, inv_dr, base_pos, d_virial, dF, sumF, f_decomp_err, param);
    // std::cout << f_decomp_err << std::endl;

    dr_new[1] = pos[0] - base_pos;
    dr_new[3] = pos[1] - base_pos;
    dr_new[4] = pos[2] - base_pos;

    show_result(dF, dr_new, base_pos, fout);

    if (base_pos[0] > up_pos[0]) fout << std::endl;
  }
}

void show_force_central(const double3* pos,
			const double3* dr,
			const double* inv_dr,
			const double3& min_pos,
			double3* dr_new,
			tensor3d& d_virial,
			double3* dF,
			double3* sumF,
			double& f_decomp_err,
			const Parameter& param,
			const std::string cur_dir) {
  auto fname = cur_dir + "/fcent.txt";
  std::ofstream fout(fname.c_str());
  fout << std::setprecision(10);
  const auto relaxed_pos = get_relaxed_pos(pos[0], pos[1], sumF[0], sumF[1]);
  const auto dr_cent = relaxed_pos - pos[1];
  
  fout << pos[0] << " " << dr[0] << std::endl;
  fout << pos[1] << " " << dr[1] << std::endl;
  fout << pos[0] << " " << dr[0] + dr[1] << std::endl;

  fout << pos[1] << " " << dr_cent << "# not minimum " << std::endl;
  fout << pos[1] << " " << min_pos - pos[1] << "# minimum stress pos" << std::endl;
  
    // compare with central force decomposition.
  fname = cur_dir + "/comp_central_force.txt";
  fout.close();
  fout.open(fname.c_str());
  fout << std::setprecision(10);
  calc_cosine_bend(pos, dr, inv_dr, relaxed_pos, d_virial, dF, sumF, f_decomp_err, param);
  dr_new[1] = pos[0] - relaxed_pos;
  dr_new[3] = pos[1] - relaxed_pos;
  dr_new[4] = pos[2] - relaxed_pos;
  fout << calc_stress_sum(dF, dr_new, 5) << " ";
  std::cout << f_decomp_err << std::endl;
  
  const double3 pos_in_line = (pos[0] + pos[2]) * 0.5;
  calc_cosine_bend(pos, dr, inv_dr, pos_in_line, d_virial, dF, sumF, f_decomp_err, param); // central force
  std::cout << "Central force decomposition check: " << dF[3].norm2() << std::endl;
  dr_new[1] = pos[0] - pos_in_line;
  dr_new[3] = pos[1] - pos_in_line;
  dr_new[4] = pos[2] - pos_in_line;
  fout << calc_stress_sum(dF, dr_new, 5) << std::endl;
}

void do_analysis(std::array<double3, 3>& pos,
		 const Parameter& param,
		 const std::string& cur_dir,
		 const std::string mode) {
  const double3 dr[] {
    pos[1] - pos[0],
    pos[2] - pos[1],
    pos[0] - pos[2]
  };
  
  const double dist[] {
    dr[0].dist2(),
    dr[1].dist2()
  };
  
  const double inv_dr[] {
    1.0 / std::sqrt(dist[0]),
    1.0 / std::sqrt(dist[1])
  };
  
  double3 dF[] {
    double3(0.0, 0.0, 0.0), // 2 -> 1
    double3(0.0, 0.0, 0.0), // 4 -> 1
    double3(0.0, 0.0, 0.0), // 3 -> 2
    double3(0.0, 0.0, 0.0), // 4 -> 2
    double3(0.0, 0.0, 0.0)  // 4 -> 3
  };
  
  double3 dr_new[] {
    double3(0.0, 0.0, 0.0), // 2 -> 1
    double3(0.0, 0.0, 0.0), // 4 -> 1
    double3(0.0, 0.0, 0.0), // 3 -> 2
    double3(0.0, 0.0, 0.0), // 4 -> 2
    double3(0.0, 0.0, 0.0)  // 4 -> 3
  };

  dr_new[0] = pos[0] - pos[1];
  dr_new[2] = pos[1] - pos[2];
  
  double3 sumF[] {
    double3(0.0, 0.0, 0.0),
    double3(0.0, 0.0, 0.0),
    double3(0.0, 0.0, 0.0)
  };

  tensor3d d_virial = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double f_decomp_err = 0.0;
  
  const auto min_pos = get_stress_sum_min_pos(dr, pos[1]);
  if (mode == modes[0]) {
    analyze_at_minimum_stress_pos(&pos[0], dr, inv_dr, min_pos, dr_new, d_virial, dF, sumF, f_decomp_err, param);    
  } else if (mode == modes[1]) {
    std::string fname = cur_dir + "/f_decomp_info.txt";
    std::ofstream fout(fname.c_str());
    analyze_along_bisect(&pos[0], dr, inv_dr, min_pos, dr_new, d_virial, dF, sumF, f_decomp_err, param, fout);
  } else if (mode == modes[2]) {
    std::string fname = cur_dir + "/stress_dist.txt";
    std::ofstream fout(fname.c_str());
    
    analyze_all_of_plane(&pos[0], dr, inv_dr, dr_new, d_virial, dF, sumF, f_decomp_err, param, fout);
    show_force_central(&pos[0], dr, inv_dr, min_pos, dr_new, d_virial, dF, sumF, f_decomp_err, param, cur_dir);
  } else {
    std::cerr << mode << ": Unknown mode.\n";
    std::exit(1);
  }
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "argv[1] == target directory name.\n";
    std::cerr << "argv[2] == execution mode name.\n";
    std::cout << modes[0] << " " << modes[1] << " " << modes[2] << std::endl;
    std::exit(1);
  }

  Parameter param(argv[1]);
  param.LoadParam();
  param.LoadCheck();

  LSParam lsparam;
  lsparam.LoadParam(argv[1]);
  
  std::array<double3, 3> pos;
  gen_init_config(pos, lsparam.bleng, lsparam.bend_angle);

  const std::string mode = std::string(argv[2]);
  
  const std::string cur_dir = std::string(argv[1]);
  do_analysis(pos, param, cur_dir, mode);
}
