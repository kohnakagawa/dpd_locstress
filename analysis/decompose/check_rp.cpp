#include "ls_tools.hpp"
#include <iostream>
#include <memory>
#include <cstdlib>
#include <sstream>
#include <numeric>
#include <limits>
#include <iomanip>
#include <random>
#include <initializer_list>
#include <algorithm>

int Parameter::sys_size;
const double offset_len = 1000.0;

void rotate_by_z(std::array<double3, 3>& pos,
		 const double angle) {
  for (auto& p : pos) {
    const auto old_pos = p;
    p.x = std::cos(angle) * old_pos.x - std::sin(angle) * old_pos.y;
    p.y = std::sin(angle) * old_pos.x + std::cos(angle) * old_pos.y;
  }
}

void gen_init_config(std::array<double3, 3>& pos,
		     const std::array<double, 2> b_len,
		     const double theta) {
  pos[0].x = pos[0].y = pos[0].z = 0.0;
  pos[1] = pos[0];
  pos[1].x += b_len[0];
  pos[2] = pos[1] + double3(b_len[1] * std::cos(theta), b_len[1] * std::sin(theta), 0.0);

  const auto rot_angle = -std::atan(pos[2].y / pos[2].x);
  rotate_by_z(pos, rot_angle);
  const double3 offset(offset_len, offset_len, 0.0);
  for (auto& p : pos) p += offset;
}

void calc_cosine_bend(const double3* r,
		      const double3* dr,
		      const double*  inv_dr,
		      const double3& base_pos,
		      double3* dF,
		      double3* sumF,
		      double& f_decomp_err,
		      double& stress_yy,
		      const Parameter& param) {
  const double	inv_dr_prod	= inv_dr[0] * inv_dr[1];
  const double	inv_dist[2]	= {inv_dr[0] * inv_dr[0],
				   inv_dr[1] * inv_dr[1]};
  const double	in_prod		= dr[0] * dr[1];
  const double	cf_b		= param.GetIntractions().cf_bend * inv_dr_prod;
  const double  cf_crs[2]       = {in_prod * inv_dist[0],
				   in_prod * inv_dist[1]};

  const double3 Ftb0(cf_b * (dr[1].x - cf_crs[0] * dr[0].x),
		     cf_b * (dr[1].y - cf_crs[0] * dr[0].y),
		     cf_b * (dr[1].z - cf_crs[0] * dr[0].z));
  const double3 Ftb1(cf_b * (dr[0].x - cf_crs[1] * dr[1].x),
		     cf_b * (dr[0].y - cf_crs[1] * dr[1].y),
		     cf_b * (dr[0].z - cf_crs[1] * dr[1].z));
  
  const double3 Ftb_sum = Ftb0 - Ftb1;
  
  decompose_3n(r[0], r[1], r[2], base_pos, -Ftb0, Ftb_sum, Ftb1, f_decomp_err, dF);

  stress_yy += dr[0].y * Ftb0.y + dr[1].y * Ftb1.y;
  
  sumF[0] = -Ftb0;
  sumF[1] = Ftb_sum;
  sumF[2] = Ftb1;
}

void dist_pair_force_stress(const double3& rj,
			    const double3& drji,
			    const double3& dFji,
			    double grid_len,
			    const double box_len,
			    std::vector<double>& buf_ls) {
  int grid_num = static_cast<int>(box_len / grid_len);
  grid_len = box_len / grid_num;
  
  const auto ri = rj + drji;
  const int i_grid = static_cast<int>(std::floor(ri.x / grid_len));
  const int j_grid = static_cast<int>(std::floor(rj.x / grid_len));

  const double stress_yy = drji.y * dFji.y;
  
  if (j_grid == i_grid) {
    buf_ls[i_grid] += stress_yy;
  } else {
    auto diff_grid = i_grid - j_grid;
    
    std::vector<double> ls_lambda;
    ls_lambda.push_back(0.0);
    if (diff_grid > 0) {
      diff_grid++;
      for (int i = 1; i < diff_grid; i++) {
	const auto wall_pos = (i + j_grid) * grid_len;
	ls_lambda.push_back((wall_pos - rj.x) / drji.x);
      }
    } else {
      diff_grid++;
      for (int i = diff_grid; i <= 0; i++) {
	const auto wall_pos = (i + j_grid) * grid_len;
	ls_lambda.push_back((wall_pos - rj.x) / drji.x);
      }
    }
    ls_lambda.push_back(1.0);

    std::sort(ls_lambda.begin(), ls_lambda.end());
    
    const int num_spreaded_cell = ls_lambda.size();
    for (int i = 1; i < num_spreaded_cell; i++) {
      assert(ls_lambda[i - 1] >= 0.0);
      assert(ls_lambda[i - 1] <= 1.0);
      
      double base_pos = rj.x + drji.x * 0.5 * (ls_lambda[i - 1] + ls_lambda[i]);
      const int base_grid = static_cast<int>(std::floor(base_pos / grid_len));
      const auto d_lambda = ls_lambda[i] - ls_lambda[i - 1];
      buf_ls[base_grid] += stress_yy * d_lambda;
    }
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

  LSParam lsparam;
  lsparam.LoadParam(argv[1]);
  
  std::array<double3, 3> pos, force;
  gen_init_config(pos, lsparam.bleng, lsparam.bend_angle);
  force.fill(double3(0.0, 0.0, 0.0));

  const double3 dr[] {
    pos[1] - pos[0],
    pos[2] - pos[1],
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

  // const double lambda = 100.0;
  // double3 base_pos = pos[1];
  // base_pos.y += lambda;
  double3 base_pos = (pos[0] + pos[2]) * 0.5;
  
  double f_decomp_err = 0.0, virial = 0.0;
  calc_cosine_bend(&pos[0], &dr[0], &inv_dr[0], base_pos, dF, &force[0], f_decomp_err, virial, param);
  
  const auto relaxed_pos = pos[1] - force[1] * (dr[0] * dr[0] / (force[1] * dr[0]));
  
  std::cout << "particle positions (" << pos[0] << ") (" << pos[1] << ") (" << pos[2] << ")" << std::endl;
  std::cout << "relaxed_pos " << relaxed_pos << std::endl;
  
  dr_new[1] = pos[0] - base_pos;
  dr_new[3] = pos[1] - base_pos;
  dr_new[4] = pos[2] - base_pos;
  
  double grid_len = 0.1, box_len = 1000.0 + offset_len;
  const int grid_num = static_cast<int>(box_len / grid_len);
  std::vector<double> buf_ls(grid_num, 0.0), buf_ls_rp(grid_num, 0.0);
  dist_pair_force_stress(pos[1], dr_new[0], dF[0], grid_len, box_len, buf_ls);
  dist_pair_force_stress(base_pos, dr_new[1], dF[1], grid_len, box_len, buf_ls);
  dist_pair_force_stress(pos[2], dr_new[2], dF[2], grid_len, box_len, buf_ls);
  dist_pair_force_stress(base_pos, dr_new[3], dF[3], grid_len, box_len, buf_ls);
  dist_pair_force_stress(base_pos, dr_new[4], dF[4], grid_len, box_len, buf_ls);

  dr_new[1] = pos[0] - relaxed_pos;
  dr_new[3] = pos[1] - relaxed_pos;
  dr_new[4] = pos[2] - relaxed_pos;
  
  dist_pair_force_stress(relaxed_pos, dr_new[1], force[0], grid_len, box_len, buf_ls_rp);
  dist_pair_force_stress(relaxed_pos, dr_new[3], force[1], grid_len, box_len, buf_ls_rp);
  dist_pair_force_stress(relaxed_pos, dr_new[4], force[2], grid_len, box_len, buf_ls_rp);

  std::cout << "total virial " << virial << std::endl;
  std::cout << "total virial buf_ls " << std::accumulate(buf_ls.cbegin(), buf_ls.cend(), 0.0) << std::endl;
  std::cout << "total virial buf_ls_rp " << std::accumulate(buf_ls_rp.cbegin(), buf_ls_rp.cend(), 0.0) << std::endl;
  
  for (size_t i = 0; i < buf_ls.size(); i++) {
    if ((buf_ls[i] != 0.0) || (buf_ls_rp[i] != 0.0))
      std::cout << (i + 0.5) * grid_len << " " << buf_ls[i] << " " << buf_ls_rp[i] << std::endl;
  }
}
