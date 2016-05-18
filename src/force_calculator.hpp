#pragma once

#include "parameter.hpp"
#include "allocate.hpp"
#include "rng.hpp"
#include <numeric>
#include <omp.h>

#define INLINE __attribute__((always_inline))

class F_calculator{
  Interactions itrs;
  BindInfo binfo;
  std::vector<double3> buf_vir;
  std::vector<double>  buf_lap_pot;
  
  int grid_numbtw = -1, numb_band = -1, p_num_band[3] = {-1, -1, -1};

  INLINE void StoreBondForce(const double3& __restrict dr,
			     const double& __restrict inv_dr,
			     double3& d_virial,
			     double& lap_conf,
			     double3* __restrict F) {
    const double cf_bond = itrs.cf_spring * (inv_dr - Parameter::i_bleng);
    
    const double3 Fbond(cf_bond * dr.x, cf_bond * dr.y, cf_bond * dr.z);
    
    d_virial.x += dr.x * Fbond.x;
    d_virial.y += dr.y * Fbond.y;
    d_virial.z += dr.z * Fbond.z;

    lap_conf += itrs.cf_spring * (6.0 * Parameter::i_bleng - 4.0 * inv_dr);
    
    F[0] -= Fbond;
    F[1] += Fbond;
  }

  INLINE void StoreBendForce(const double3* __restrict dr,
			     const double * __restrict inv_dr,
			     const double * __restrict dist,
			     double3& d_virial,
			     double& lap_conf,
			     double3* __restrict F) {
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
    
    d_virial.x += dr[0].x * Ftb0.x + dr[1].x * Ftb1.x;
    d_virial.y += dr[0].y * Ftb0.y + dr[1].y * Ftb1.y;
    d_virial.z += dr[0].z * Ftb0.z + dr[1].z * Ftb1.z;

    lap_conf += 2.0 * cf_b * inv_dist[0] * inv_dist[1] * (in_prod * (in_prod + 2.0 * (dist[0] + dist[1])) + dist[0] * dist[1]);    
    
    F[0] -= Ftb0;
    F[1] += Ftb0 - Ftb1;
    F[2] += Ftb1;
  }

  INLINE void ClearForce(double3* force) {
    for (int i = 0; i < Parameter::sys_size; i++) force[i].clear();
    buf_lap_pot.assign(buf_lap_pot.size(), 0.0);
    buf_vir.assign(buf_vir.size(), double3(0.0));
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
		    double3& __restrict d_virial,
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
  
public:
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

  double3 DumpVirial() const {
    return std::accumulate(buf_vir.cbegin(), buf_vir.cend(), double3(0.0, 0.0, 0.0));
  }
};

#undef INLINE
