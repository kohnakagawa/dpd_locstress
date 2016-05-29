#pragma once

#include "parameter.hpp"
#include "allocate.hpp"
#include "rng.hpp"
#include "lapack_wrapper.hpp"
#include <numeric>
#include <omp.h>
#include <algorithm>
#include <iomanip>

#ifdef __INTEL_COMPILER
#define CALL_DGELSD call_mkl_dgelsd
#else
#define CALL_DGELSD call_f77_dgelsd
#endif

#define INLINE __attribute__((always_inline))

class F_calculator {
public:
  Interactions itrs;
  BindInfo binfo;
  std::vector<tensor3d> buf_vir;
  std::vector<double>  buf_lap_pot;

  static constexpr int ls_type_pot = NUM_LS_TYPE - 1;
  
  std::array<std::vector<std::vector<tensor3d> >, ls_type_pot> buf_lstress;
  double f_decomp_err = 0.0, accum_stress = 0.0;
  tensor3d vir_sum;
  
  int grid_numbtw = -1, numb_band = -1, p_num_band[3] = {-1, -1, -1};

  static INLINE std::array<int, 3> GetLSGrid(const double3& r,
					     const Parameter& param) {
    const std::array<int, 3> grid_id = {
      static_cast<int>(std::floor(r.x * param.i_ls_grid_.x)),
      static_cast<int>(std::floor(r.y * param.i_ls_grid_.y)),
      static_cast<int>(std::floor(r.z * param.i_ls_grid_.z))
    };
    return grid_id;
  }

  static INLINE int GetLSGrid1d(std::array<int, 3> grid_id,
				const Parameter& param) {
    if (grid_id[0] < 0 || grid_id[0] >= param.ls_grid_num_[0])
      grid_id[0] -= static_cast<int>(grid_id[0] * param.ls_grid_num_[0] / std::abs(grid_id[0]));

    if (grid_id[1] < 0 || grid_id[1] >= param.ls_grid_num_[1])
      grid_id[1] -= static_cast<int>(grid_id[1] * param.ls_grid_num_[1] / std::abs(grid_id[1]));

    if (grid_id[2] < 0 || grid_id[2] >= param.ls_grid_num_[2])
      grid_id[2] -= static_cast<int>(grid_id[2] * param.ls_grid_num_[2] / std::abs(grid_id[2]));

    const auto ret = grid_id[0] + param.ls_grid_num_[0] * (grid_id[1] + grid_id[2] * param.ls_grid_num_[1]);
#ifdef DEBUG
    CHECK_EQUATION(grid_id[0] < param.ls_grid_num_[0], grid_id[0]);
    CHECK_EQUATION(grid_id[1] < param.ls_grid_num_[1], grid_id[1]);
    CHECK_EQUATION(grid_id[2] < param.ls_grid_num_[2], grid_id[2]);
    CHECK_EQUATION(ret < param.ls_grid_num_[0] * param.ls_grid_num_[1] * param.ls_grid_num_[2], ret);
#endif
    return ret;
  }

  // NOTE: this function is applied for each axis.
  INLINE void CalcLSLambdaForEachAxis(std::vector<double>& lambda,
				      const int j_grid,
				      int diff_grid,
				      const double rj,
				      const double drji,
				      const int axis,
				      const Parameter& param) const {
    if (diff_grid > 0) {
      diff_grid++;
      for (int i = 1; i < diff_grid; i++) {
	const double wall_pos = (i + j_grid) * param.ls_grid_[axis];
	lambda.push_back((wall_pos - rj) / drji);
      }
    } else if (diff_grid < 0) {
      diff_grid++;
      for (int i = diff_grid; i <= 0; i++) {
	const double wall_pos = (i + j_grid) * param.ls_grid_[axis];
	lambda.push_back((wall_pos - rj) / drji);
      }
    }
  }

  INLINE std::vector<double> CalcLSLambda(const std::array<int, 3>& j_grid,
					  const std::array<int, 3>& diff_grid,
					  const double3 rj,
					  const double3 drji,
					  const Parameter& param) {
    std::vector<double> lambda;
    lambda.push_back(0.0);
    CalcLSLambdaForEachAxis(lambda, j_grid[0], diff_grid[0], rj[0], drji[0], 0, param);
    CalcLSLambdaForEachAxis(lambda, j_grid[1], diff_grid[1], rj[1], drji[1], 1, param);
    CalcLSLambdaForEachAxis(lambda, j_grid[2], diff_grid[2], rj[2], drji[2], 2, param);
    lambda.push_back(1.0);
    std::sort(lambda.begin(), lambda.end());
    return lambda;
  }

  static INLINE void ApplyPBC(double3& r, const Parameter& param) {
    r.x -= std::floor(r.x * param.iL.x) * param.L.x;
    r.y -= std::floor(r.y * param.iL.y) * param.L.y;
    r.z -= std::floor(r.z * param.iL.z) * param.L.z;
  }
  
  // NOTE: this function is thread safe.
  void DistPairForceStress(const double3& rj,
			   const double3& drji,
			   const double3& dFji,
			   const Parameter& param,
			   const int stress_type) {
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
	  buf_lstress[stress_type][tid][j_grid_1d][i][j] += stress[i][j];
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
	    buf_lstress[stress_type][tid][base_grid][j][k] += stress[j][k] * d_lambda;
      }
    }
  }

  // NOTE: Central force decomposition
  void Decompose3NCfd(const double3& r1, const double3& r2, const double3& r3,
		      const double3& F1, const double3& F2, const double3& F3,
		      const Parameter& param) {
    auto dr13 = r3 - r1; MinImage(dr13, param);
    auto base_pos = r1 + 0.5 * dr13; ApplyPBC(base_pos, param);
    Decompose3N(r1, r2, r3, base_pos, F1, F2, F3, param);
  }
  
  void Decompose3N(const double3& r1, const double3& r2, const double3& r3,
		   const double3& F1, const double3& F2, const double3& F3,
		   const Parameter& param) {
    auto dr21 = r1 - r2; MinImage(dr21, param);
    auto cross_point = r2 + F2 * (dr21 * dr21 / (F2 * dr21));

    auto dr42 = r2 - cross_point;
    auto dr41 = dr42 + dr21;
    auto dr23 = r3 - r2; MinImage(dr23, param);
    auto dr43 = dr42 + dr23;

    ApplyPBC(cross_point, param);
    
    DistPairForceStress(cross_point, dr41, F1, param, ANGLE);
    DistPairForceStress(cross_point, dr42, F2, param, ANGLE);
    DistPairForceStress(cross_point, dr43, F3, param, ANGLE);
  }

  // General case
  // NOTE: this function is NOT thread safe.
  void Decompose3N(const double3& r1, const double3& r2, const double3& r3, const double3& r4, // in-plane (ri-rj ^ ri-rk) position
		   const double3& F1, const double3& F2, const double3& F3,
		   const Parameter& param) {
    auto dr21 = r1 - r2; MinImage(dr21, param);
    auto dr41 = r1 - r4; MinImage(dr41, param);
    auto dr32 = r2 - r3; MinImage(dr32, param);
    auto dr42 = r2 - r4; MinImage(dr42, param);
    auto dr43 = r3 - r4; MinImage(dr43, param);

    constexpr int nRows = 12, nCols = 5, nRhs = 1;

    std::array<double, nRows * nCols> D; D.fill(0.0);
    std::array<double, nRows> b;

    D[nRows * 0 + 0] =  dr21.x; D[nRows * 1 + 0] = dr41.x;
    D[nRows * 0 + 1] =  dr21.y; D[nRows * 1 + 1] = dr41.y;
    D[nRows * 0 + 2] =  dr21.z; D[nRows * 1 + 2] = dr41.z;
    b[0] = F1.x; b[1] = F1.y; b[2] = F1.z;
    
    D[nRows * 0 + 3] = -dr21.x; D[nRows * 2 + 3] = dr32.x; D[nRows * 3 + 3] = dr42.x;
    D[nRows * 0 + 4] = -dr21.y; D[nRows * 2 + 4] = dr32.y; D[nRows * 3 + 4] = dr42.y;
    D[nRows * 0 + 5] = -dr21.z; D[nRows * 2 + 5] = dr32.z; D[nRows * 3 + 5] = dr42.z;
    b[3] = F2.x; b[4] = F2.y; b[5] = F2.z;
    
    D[nRows * 2 + 6] = -dr32.x; D[nRows * 4 + 6] = dr43.x;
    D[nRows * 2 + 7] = -dr32.y; D[nRows * 4 + 7] = dr43.y;
    D[nRows * 2 + 8] = -dr32.z; D[nRows * 4 + 8] = dr43.z;
    b[6] = F3.x; b[7] = F3.y; b[8] = F3.z;
    
    D[nRows * 1 + 9]  = -dr41.x; D[nRows * 3 + 9]  = -dr42.x; D[nRows * 4 + 9]  = -dr43.x;
    D[nRows * 1 + 10] = -dr41.y; D[nRows * 3 + 10] = -dr42.y; D[nRows * 4 + 10] = -dr43.y;
    D[nRows * 1 + 11] = -dr41.z; D[nRows * 3 + 11] = -dr42.z; D[nRows * 4 + 11] = -dr43.z;
    b[9] = 0.0;  b[10] = 0.0; b[11] = 0.0;

    CALL_DGELSD<nRows, nCols, nRhs>(D, b);

    const double3 dF21(b[0] * dr21.x, b[0] * dr21.y, b[0] * dr21.z);
    const double3 dF41(b[1] * dr41.x, b[1] * dr41.y, b[1] * dr41.z);
    const double3 dF32(b[2] * dr32.x, b[2] * dr32.y, b[2] * dr32.z);
    const double3 dF42(b[3] * dr42.x, b[3] * dr42.y, b[3] * dr42.z);
    const double3 dF43(b[4] * dr43.x, b[4] * dr43.y, b[4] * dr43.z);

    // check err
    const double3 F1_d = dF21 + dF41;
    const double3 F2_d = dF32 + dF42 - dF21;
    const double3 F3_d = dF43 - dF32;
    const double3 F4_d = -dF41 - dF42 - dF43;
    
    f_decomp_err += (F1_d - F1) * (F1_d - F1);
    f_decomp_err += (F2_d - F2) * (F2_d - F2);
    f_decomp_err += (F3_d - F3) * (F3_d - F3);
    f_decomp_err += F4_d * F4_d;

    accum_stress += dF21.norm2() * dr21.norm2()
      + dF41.norm2() * dr41.norm2()
      + dF32.norm2() * dr32.norm2()
      + dF42.norm2() * dr42.norm2()
      + dF43.norm2() * dr43.norm2();
    
    DistPairForceStress(r2, dr21, dF21, param, ANGLE);
    DistPairForceStress(r4, dr41, dF41, param, ANGLE);
    DistPairForceStress(r3, dr32, dF32, param, ANGLE);
    DistPairForceStress(r4, dr42, dF42, param, ANGLE);
    DistPairForceStress(r4, dr43, dF43, param, ANGLE);
  }

  INLINE void StoreBondForce(const double3& __restrict r,
			     const double3& __restrict dr,
			     const double&  __restrict inv_dr,
			     tensor3d&      __restrict d_virial,
			     double&        __restrict lap_conf,
			     double3*       __restrict F,
			     const Parameter& __restrict param) {
    const double cf_bond = itrs.cf_spring * (inv_dr - Parameter::i_bleng);
    
    const double3 Fbond(cf_bond * dr.x, cf_bond * dr.y, cf_bond * dr.z);
    
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) d_virial[i][j] += dr[i] * Fbond[j];

    lap_conf += itrs.cf_spring * (6.0 * Parameter::i_bleng - 4.0 * inv_dr);
    
    F[0] -= Fbond;
    F[1] += Fbond;

#ifdef CALC_LOC_STRESS
    DistPairForceStress(r, dr, Fbond, param, BOND);
#endif
  }

  INLINE void StoreBendForce(const double3* __restrict r,
			     const double3* __restrict dr,
			     const double*  __restrict inv_dr,
			     const double*  __restrict dist,
			     tensor3d&      __restrict d_virial,
			     double&        __restrict lap_conf,
			     double3*       __restrict F,
			     const Parameter& __restrict param) {
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
    
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
	d_virial[i][j] += dr[0][i] * Ftb0[j] + dr[1][i] * Ftb1[j];	
      }

    lap_conf += 2.0 * cf_b * (in_prod * (2.0 * (inv_dist[0] + inv_dist[1]) + in_prod * inv_dr_prod * inv_dr_prod) + 1.0);
    
    const double3 Ftb_sum = Ftb0 - Ftb1;
    
    F[0] -= Ftb0;
    F[1] += Ftb_sum;
    F[2] += Ftb1;

#ifdef CALC_LOC_STRESS
#ifdef CENTRAL_FORCE
    Decompose3NCfd(r[0], r[1], r[2], -Ftb0, Ftb_sum, Ftb1, param);
#elif defined RELAXED_BASE_POS
    Decompose3N(r[0], r[1], r[2], -Ftb0, Ftb_sum, Ftb1, param);
#else // general case
    const double3 bi_vec = Ftb_sum / Ftb_sum.norm2();
    double3 base_pos = r[1] + bi_vec * param.ls_lambda;
    ApplyPBC(base_pos, param);
    Decompose3N(r[0], r[1], r[2], base_pos, -Ftb0, Ftb_sum, Ftb1, param);
#endif // end of CENTRAL_FORCE
#endif // end of CALC_LOC_STRESS
  }

  INLINE void ClearForce(double3* force) {
    for (int i = 0; i < Parameter::sys_size; i++) force[i].clear();
    buf_lap_pot.assign(buf_lap_pot.size(), 0.0);
    buf_vir.assign(buf_vir.size(), tensor3d{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}});

    for (int i = 0; i < ls_type_pot; i++) {
      const int num_lstress_grid = buf_lstress[i].size();
      for (int j = 0; j < num_lstress_grid; j++) {
	buf_lstress[i][j].assign(buf_lstress[i][j].size(), tensor3d{{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}});
      }
    }
    f_decomp_err = 0.0;
    accum_stress = 0.0;
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
    
    CHECK_EQUATION(param.grid_numb[2] % numb_band == 0, param.grid_numb[2]);
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

  void Collision_pow_one(double3& __restrict vi,
			 double3& __restrict vj,
			 const double3& __restrict dr,
			 const double cf_g,
			 const double nrml,
			 const double inv_dr,
			 const Parameter& __restrict param) {
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

  void Collision_pow_half(double3& __restrict vi,
			  double3& __restrict vj,
			  const double3& __restrict dr,
			  const double cf_g,
			  const double nrml,
			  const double inv_dr,
			  const double dr_norm, 
			  const Parameter& __restrict param) {
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

  void AddCollisionHalf(const double3*   __restrict pr,
			double3*         __restrict pv,
			const par_prop*  __restrict prop,
			const CellList&  __restrict clist,
			const Parameter& __restrict param,
			RNG&             __restrict rng,
			const int                   call_num);

  //calculate force in amphiphile molecule
  void CalcForceInAmp(const double3* __restrict pr,
		      double3*       __restrict force,
		      const ChemInfo& __restrict cheminfo,
		      const Parameter& __restrict param);
  
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
    
    StoreBondForce(temp_pos[0], dr[0], inv_dr[0], d_virial, lap_conf, &Fbb[0], param);

    for (int unit = 2; unit < bond_n; unit++) {
      const int load_dest = elem_idx[beg_idx + unit];
      temp_pos[unit] = pr[load_dest];
      dr[unit - 1] = temp_pos[unit] - temp_pos[unit - 1];
      MinImage(dr[unit - 1], param);
      dist2[unit - 1] = dr[unit - 1].dist2();
      inv_dr[unit - 1] = 1.0 / std::sqrt(dist2[unit - 1]);
      
      StoreBondForce(temp_pos[unit - 1], dr[unit - 1], inv_dr[unit - 1], d_virial, lap_conf, &Fbb[unit - 1], param);
      StoreBendForce(&temp_pos[unit - 2], &dr[unit - 2], &inv_dr[unit - 2], dist2, d_virial, lap_conf, &Fbb[unit - 2], param);
    }

    for (int unit = 0; unit < bond_n; unit++) {
      const int str_dest = elem_idx[beg_idx + unit];
      force[str_dest] += Fbb[unit];
    }
  }
  
  explicit F_calculator(const Parameter& param) {
    itrs  = param.GetIntractions();
    binfo = param.GetBindInform();
    const int th_numb   = omp_get_max_threads();
    buf_vir.resize(th_numb, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    buf_lap_pot.resize(th_numb, 0.0);
    
    const auto all_ls_grid = param.ls_grid_num_[0] * param.ls_grid_num_[1] * param.ls_grid_num_[2];
    for (int i = 0; i < ls_type_pot; i++) {
      buf_lstress[i].resize(th_numb);
      for (int j = 0; j < th_numb; j++) {
	buf_lstress[i][j].resize(all_ls_grid, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
      }
    }
    vir_sum.fill({0.0, 0.0, 0.0});
    DevideCell(param);
  }
  
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

  const tensor3d& DumpVirial() {
    vir_sum.fill({0.0, 0.0, 0.0});
    const int th_numb = omp_get_max_threads();
    for (int tid = 0; tid < th_numb; tid++)
      for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) {
	  vir_sum[j][k] += buf_vir[tid][j][k];
	}
    return vir_sum;
  }

  tensor3d DumpVirialError() const {
    tensor3d vir_sum_ls = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    const int thnum = buf_lstress[0].size();
    const int grid_num = buf_lstress[0][0].size();
    
    for (int t = 0; t < ls_type_pot; t++)
      for (int i = 0; i < thnum; i++)
	for (int g = 0; g < grid_num; g++)
	  for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) {
	      vir_sum_ls[j][k] += buf_lstress[t][i][g][j][k];
	  }
    return vir_sum_ls - vir_sum;
  }

  double DumpFdecompError() const {
    return std::sqrt(f_decomp_err);
  }

  double DumpStressAccumu() const {
    return accum_stress;
  }

  const std::vector<std::vector<tensor3d> >& DumpCurLocStress(const int type) const {
    return buf_lstress[type];
  }
};

#undef INLINE
