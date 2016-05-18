#pragma once

#include "parameter.hpp"
#include "allocate.hpp"
#include "rng.hpp"
#include "Eigen/Dense"
#include <numeric>
#include <omp.h>
#include <algorithm>

#define INLINE __attribute__((always_inline))

class F_calculator{
  Interactions itrs;
  BindInfo binfo;
  std::vector<tensor3d> buf_vir;
  std::vector<double>  buf_lap_pot;
  std::vector<std::vector<tensor3d> > buf_lstress;
  
  int grid_numbtw = -1, numb_band = -1, p_num_band[3] = {-1, -1, -1};

  inline std::array<int, 3> GetLSGrid(const double3& r,
				      const Parameter& param) const {
    const std::array<int, 3> grid_id = {
      static_cast<int>(r.x * param.i_ls_grid_.x),
      static_cast<int>(r.y * param.i_ls_grid_.y),
      static_cast<int>(r.z * param.i_ls_grid_.z)
    };
    return grid_id;
  }

  inline int GetLSGrid1d(const std::array<int, 3>& grid_id,
			 const Parameter& param) const {
    const auto ret = grid_id[0] + param.ls_grid_num_[0] * (grid_id[1] + grid_id[2] * param.ls_grid_num_[1]);
#ifdef DEBUG
    CHECK_EQUATION(ret < param.ls_grid_num_[0] * param.ls_grid_num_[1] * param.ls_grid_num_[2], ret);
#endif
    return ret;
  }

  // NOTE: this function is applied for each axis.
  inline void CalcLSLambdaForEachAxis(std::vector<double>& lambda,
				      const int j_grid,
				      const int diff_grid,
				      const double rj_pos,
				      const double drji,
				      const int axis,
				      const Parameter& param) const {
    for (int i = 0; i < diff_grid; i++) {
      const double wall_pos = (i + j_grid) * param.ls_grid_[axis];
      lambda.push_back((wall_pos - rj_pos) / drji);
    }
  }

  inline std::vector<double> CalcLSLambda(const std::array<int, 3>& j_grid,
					  const std::array<int, 3>& diff_grid,
					  const double3 rj_pos,
					  const double3 drji,
					  const Parameter& param) {
    std::vector<double> lambda = {0.0};
    CalcLSLambdaForEachAxis(lambda, j_grid[0], diff_grid[0], rj_pos[0], drji[0], 0, param);
    CalcLSLambdaForEachAxis(lambda, j_grid[1], diff_grid[1], rj_pos[1], drji[1], 1, param);
    CalcLSLambdaForEachAxis(lambda, j_grid[2], diff_grid[2], rj_pos[2], drji[2], 2, param);
    lambda.push_back(1.0);
    std::sort(lambda.begin(), lambda.end());
    return lambda;
  }

  void ApplyPBC(double3& r, const Parameter& param) {
    r.x -= std::floor(r.x * param.iL.x) * param.L.x;
    r.y -= std::floor(r.y * param.iL.y) * param.L.y;
    r.z -= std::floor(r.z * param.iL.z) * param.L.z;
  }
  
  void DistPairForceStress(const double3& rj,
			   const double3& drji,
			   const double3& dFji,
			   const Parameter& param) {
    const int tid = omp_get_thread_num();

    const auto i_grid = GetLSGrid(rj + drji, param);
    const auto j_grid = GetLSGrid(rj, param);

    tensor3d stress;
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
	stress[i][j] = drji[i] * dFji[j];

    if (j_grid == i_grid) {
      // simply add local stress
      const auto j_grid_1d = GetLSGrid1d(j_grid, param);
      for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
	  buf_lstress[tid][j_grid_1d][i][j] += stress[i][j];
    } else {
      // spread local stress
      const std::array<int, 3> diff_grid = {
	i_grid[0] - j_grid[0],
	i_grid[1] - j_grid[1],
	i_grid[2] - j_grid[2]
      };
      
      const auto lambda = CalcLSLambda(j_grid, diff_grid, rj, drji, param);
      const int num_spreaded_cell = lambda.size();
      for (int i = 1; i < num_spreaded_cell; i++) {
	auto base_pos = rj + drji * 0.5 * (lambda[i - 1] + lambda[i]);
	ApplyPBC(base_pos, param);
	const auto base_grid = GetLSGrid1d(GetLSGrid(base_pos, param), param);
	const auto d_lambda = lambda[i] - lambda[i - 1];
	for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++)
	    buf_lstress[tid][base_grid][j][k] += stress[j][k] * d_lambda;
      }
    }
  }

  void Decompose3N(const double3& rl, // in-plane (ri-rj ^ ri-rk) position
		   const double3& ri, const double3& rj, const double3& rk,
		   const double3& Fi, const double3& Fj, const double3& Fk,
		   const Parameter& param) {
    // Decompose forces
    auto dr21 = rj - ri; MinImage(dr21, param); const auto dr21_n = dr21 / dr21.norm2();
    auto dr41 = rl - ri; MinImage(dr41, param); const auto dr41_n = dr41 / dr41.norm2();
    auto dr32 = rk - rj; MinImage(dr32, param); const auto dr32_n = dr32 / dr32.norm2();
    auto dr42 = rl - rj; MinImage(dr42, param); const auto dr42_n = dr42 / dr42.norm2();
    auto dr43 = rl - rk; MinImage(dr43, param); const auto dr43_n = dr43 / dr43.norm2();
    
    Eigen::MatrixXd D(12, 5);
    D << dr21_n.x, dr41_n.x, 0.0, 0.0, 0.0,
         dr21_n.y, dr41_n.y, 0.0, 0.0, 0.0,
         dr21_n.z, dr41_n.z, 0.0, 0.0, 0.0,
        -dr21_n.x, -dr32_n.x, dr42_n.x, 0.0, 0.0,
        -dr21_n.y, -dr32_n.y, dr42_n.y, 0.0, 0.0,
        -dr21_n.z, -dr32_n.z, dr42_n.z, 0.0, 0.0,
         0.0, 0.0, dr32_n.x, 0.0, -dr43_n.x,
         0.0, 0.0, dr32_n.y, 0.0, -dr43_n.y,
         0.0, 0.0, dr32_n.z, 0.0, -dr43_n.z,
         0.0, -dr41_n.x, 0.0, -dr42_n.x, -dr43_n.x,
         0.0, -dr41_n.y, 0.0, -dr42_n.y, -dr43_n.y,
         0.0, -dr41_n.z, 0.0, -dr42_n.z, -dr43_n.z;
    Eigen::VectorXd F(12);
    F << Fi.x, Fi.y, Fi.z, Fj.x, Fj.y, Fj.z, Fk.x, Fk.y, Fk.z, 0.0, 0.0, 0.0;
    Eigen::VectorXd phi = D.colPivHouseholderQr().solve(F);
    
    // Distribute local stress
    const double3 dF21(phi[0] * dr21_n.x, phi[0] * dr21_n.y, phi[0] * dr21_n.z);
    const double3 dF41(phi[1] * dr41_n.x, phi[1] * dr41_n.y, phi[1] * dr41_n.z);
    const double3 dF32(phi[2] * dr32_n.x, phi[2] * dr32_n.y, phi[2] * dr32_n.z);
    const double3 dF42(phi[3] * dr42_n.x, phi[3] * dr42_n.y, phi[3] * dr42_n.z);
    const double3 dF43(phi[4] * dr43_n.x, phi[4] * dr43_n.y, phi[4] * dr43_n.z);
    
    DistPairForceStress(ri, dr21, dF21, param);
    DistPairForceStress(ri, dr41, dF41, param);
    DistPairForceStress(rj, dr32, dF32, param);
    DistPairForceStress(rj, dr42, dF42, param);
    DistPairForceStress(rk, dr43, dF43, param);    
  }

  INLINE void StoreBondForce(const double3& __restrict dr,
			     const double& __restrict inv_dr,
			     tensor3d& d_virial,
			     double& lap_conf,
			     double3* __restrict F) {
    const double cf_bond = itrs.cf_spring * (inv_dr - Parameter::i_bleng);
    
    const double3 Fbond(cf_bond * dr.x, cf_bond * dr.y, cf_bond * dr.z);
    
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) d_virial[i][j] += dr[i] * Fbond[j];
    
    lap_conf += itrs.cf_spring * (6.0 * Parameter::i_bleng - 4.0 * inv_dr);
    
    F[0] -= Fbond;
    F[1] += Fbond;
  }

  INLINE void StoreBendForce(const double3* __restrict dr,
			     const double*  __restrict inv_dr,
			     const double*  __restrict dist,
			     tensor3d&       __restrict d_virial,
			     double&        __restrict lap_conf,
			     double3*       __restrict F) {
    const double	inv_dr_prod = inv_dr[0] * inv_dr[1];
    const double	inv_dist[2] = {inv_dr[0] * inv_dr[0],
				       inv_dr[1] * inv_dr[1]};
    const double	in_prod	    = dr[0] * dr[1];
    const double	cf_b	    = itrs.cf_bend * inv_dr_prod;
    const double        cf_crs[2]   = {in_prod * inv_dist[0],
				       in_prod * inv_dist[1]};
    
    const double3 Ftb0(cf_b * (dr[1].x - cf_crs[0] * dr[0].x),
		       cf_b * (dr[1].y - cf_crs[0] * dr[0].y),
		       cf_b * (dr[1].z - cf_crs[0] * dr[0].z));
    const double3 Ftb1(cf_b * (dr[0].x - cf_crs[1] * dr[1].x),
		       cf_b * (dr[0].y - cf_crs[1] * dr[1].y),
		       cf_b * (dr[0].z - cf_crs[1] * dr[1].z));
    
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
				  d_virial[i][j] += dr[0][i] * Ftb0[j] + dr[1][i] * Ftb1[j];

    lap_conf += 2.0 * cf_b * inv_dist[0] * inv_dist[1] * (in_prod * (in_prod + 2.0 * (dist[0] + dist[1])) + dist[0] * dist[1]);    
    
    F[0] -= Ftb0;
    F[1] += Ftb0 - Ftb1;
    F[2] += Ftb1;
  }

  INLINE void ClearForce(double3* force) {
    for (int i = 0; i < Parameter::sys_size; i++) force[i].clear();
    buf_lap_pot.assign(buf_lap_pot.size(), 0.0);
    buf_vir.assign(buf_vir.size(), {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    for (size_t i = 0; i < buf_lstress.size(); i++)
      buf_lstress[i].assign(buf_lstress[i].size(), {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  }

  void CheckForce(double3* force) {
    for (int i = 0; i < Parameter::sys_size; i++) force[i].isfinite3();
  }

  void DevideCell(const Parameter& param) {
    const int th_numb   = omp_get_max_threads();
    grid_numbtw = param.grid_numb[0] * param.grid_numb[1];
    numb_band   = param.grid_numb[2] / th_numb;
  
    p_num_band[0] = 0;
    p_num_band[1] = param.grid_numb[2] / th_numb / 2;
    p_num_band[2] = param.grid_numb[2] / th_numb - p_num_band[1];

    assert(param.grid_numb[2] % numb_band == 0);
  }

  void RegistNearIdx(const int& itr_id,
		     const int& tar_id,
		     const float& dist,
		     ChemInfo& cheminfo);

  void CountNearWater(const int& itr_id,
		      const int& tar_id,		      
		      const par_prop& i_prp,
		      const par_prop& t_prp,
		      float dist,
		      ChemInfo& cheminfo);

  void CalcForceHalf(const double3*  __restrict pr,
		     const par_prop* __restrict prop,
		     double3*        __restrict force,
		     const CellList&            clist,
		     const Parameter&           param,
		     const int                  call_num);

  void Collision_pow_one(double3& vi,
			 double3& vj,
			 const double3& dr,
			 const double cf_g,
			 const double nrml,
			 const double inv_dr,
			 const Parameter& param) {
    const double  g_invdr_dt = cf_g * (inv_dr - 1.0) * param.dt_c;
    const double  cf_numer   = inv_dr * g_invdr_dt;
    const double3 dv         = vi - vj;
    const double  dvdr       = dv * dr;
    const double  cf_A       = cf_numer * dvdr;
    const double  cf_B       = std::sqrt(cf_numer) * nrml;
    const double  col_cf     = (cf_B - cf_A) * inv_dr / (inv_dr + g_invdr_dt);
    const double3 col_impls(col_cf * dr.x, col_cf * dr.y, col_cf * dr.z);
    vi += col_impls;
    vj -= col_impls;
  }

  void Collision_pow_half(double3& vi,
			  double3& vj,
			  const double3& dr,
			  const double cf_g,
			  const double nrml,
			  const double inv_dr,
			  const double dr_norm, 
			  const Parameter& param) {
    const double w_dt			= cf_g * std::sqrt(1.0 - dr_norm) * param.dt_c;
    const double3 dv			= vi - vj;
    
    const double sqrt_w_dt		= std::sqrt(w_dt);
    const double one_p_wdt		= 1.0 + w_dt;
    const double inv_one_p_wdt		= 1.0 / one_p_wdt;
    const double inv_dr_sqrt_w_dt	= inv_dr * sqrt_w_dt;
    const double dvdr			= dv * dr;
    
    const double col_cf			= inv_one_p_wdt * (inv_dr_sqrt_w_dt * (nrml - inv_dr_sqrt_w_dt * dvdr));
    const double3 col_impls(col_cf * dr.x, col_cf * dr.y, col_cf * dr.z);
    
    vi += col_impls;
    vj -= col_impls;
  }

  void AddCollisionHalf(const double3*  __restrict pr,
			double3*        __restrict pv,
			const par_prop* __restrict prop,
			const CellList&            clist,
			const Parameter&           param,
			RNG&                       rng,
			const int                  call_num);

  //calculate force in amphiphile molecule
  void CalcForceInAmp(const double3* __restrict pr,
		      double3* __restrict force,
		      const ChemInfo& cheminfo,
		      const Parameter& param);
  
  /*TODO: This part may be optimized by using loop unrolling. */
  template<int bond_n>
  void CalcBondBend(const double3* __restrict pr,
		    double3* __restrict force,
		    tensor3d& __restrict d_virial,
		    double&  __restrict lap_conf,
		    const int beg_idx,
		    const int* elem_idx,
		    const Parameter& param) {
    double3 Fbb[bond_n], temp_pos[bond_n], dr[bond_n - 1];
    double  dist2[bond_n - 1], inv_dr[bond_n - 1];

    const int l_dst[2] = {elem_idx[beg_idx], elem_idx[beg_idx + 1]};
    temp_pos[0] = pr[l_dst[0]]; temp_pos[1] = pr[l_dst[1]];
    dr[0] = temp_pos[1] - temp_pos[0];
    MinImage(dr[0], param);
    dist2[0] = dr[0].dist2();
    inv_dr[0] = 1.0 / sqrt(dist2[0]);
    
    StoreBondForce(dr[0], inv_dr[0], d_virial, lap_conf, &Fbb[0]);

    for (int unit = 2; unit < bond_n; unit++) {
      const int load_dest = elem_idx[beg_idx + unit];
      temp_pos[unit] = pr[load_dest];
      dr[unit - 1] = temp_pos[unit] - temp_pos[unit - 1];
      MinImage(dr[unit - 1], param);
      dist2[unit - 1] = dr[unit - 1].dist2();
      inv_dr[unit - 1] = 1.0 / std::sqrt(dist2[unit - 1]);
      
      StoreBondForce(dr[unit - 1], inv_dr[unit - 1], d_virial, lap_conf, &Fbb[unit - 1]);
      StoreBendForce(&dr[unit - 2], &inv_dr[unit - 2], dist2, d_virial, lap_conf, &Fbb[unit - 2]);
    }

    for (int unit = 0; unit < bond_n; unit++) {
      const int str_dest = elem_idx[beg_idx + unit];
      force[str_dest] += Fbb[unit];
    }
  }
  
  // public:
  explicit F_calculator(const Parameter& param);
  
  void AddConservForce(const double3*  __restrict pr,
		       const par_prop* __restrict prop,
		       double3*        __restrict force,
		       const CellList& clist,
		       const Parameter& param,
		       ChemInfo& cheminfo);

  void AddDisspRandom(const double3*  __restrict pr,
		      double3*        __restrict pv,
		      const par_prop* __restrict prop,
		      const CellList&            clist,
		      const Parameter&           param,
		      RNG&                       rng);
  
  void AddBindForceCyl(const double3* __restrict pr,
		       const par_prop* __restrict prop,
		       double3* __restrict force,
		       const ChemInfo& cheminfo);

  void AddBindForceSph(const double3* __restrict pr,
		       const par_prop* __restrict prop,
		       double3* __restrict force);

  void AddPosRes(const double3* __restrict pr,
		 const double3* __restrict pr_base,
		 double3* __restrict force,
		 const bool* __restrict rest_on);

  static INLINE void MinImage(double3& a, const Parameter& param) {
    a.x -= param.L.x * std::round(a.x * param.iL.x);
    a.y -= param.L.y * std::round(a.y * param.iL.y);
    a.z -= param.L.z * std::round(a.z * param.iL.z);
  }
  
  double DumpConfigT(const double3* F) const {
    const double lap_pot = std::accumulate(buf_lap_pot.cbegin(), buf_lap_pot.cend(), 0.0);
    const double numer   = std::accumulate(F, F + Parameter::sys_size, 0.0,
					   [](const double sum, const double3& val) {return sum + val * val;});
    return numer / lap_pot;
  }

  tensor3d DumpVirial() const {
    tensor3d result = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const int th_numb = omp_get_max_threads();
    for (int tid = 0; tid < th_numb; tid++)
      for (int j = 0; j < 3; j++) 
	for (int k = 0; k < 3; k++)
	  result[j][k] += buf_vir[tid][j][k];
    return result;
  }
};

#undef INLINE
