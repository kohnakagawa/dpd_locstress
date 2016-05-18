#include "force_calculator.hpp"

template<>
void F_calculator::CalcBondBend<2>(const double3* __restrict	pr,
				   double3* __restrict		force,
				   double3& __restrict          d_virial,
				   double&  __restrict          lap_conf,
				   const int			beg_idx,
				   const int*			elem_idx,
				   const Parameter&		param) {
  const int l_dst[2] = {elem_idx[beg_idx], elem_idx[beg_idx + 1]};
  double3 dr(pr[l_dst[1]].x - pr[l_dst[0]].x,
	     pr[l_dst[1]].y - pr[l_dst[0]].y,
	     pr[l_dst[1]].z - pr[l_dst[0]].z);
  MinImage(dr, param);
  const double inv_dr = dr.invnorm2();

  const double cf_bond = itrs.cf_spring * (inv_dr - Parameter::i_bleng);
  const double3 Fbond(cf_bond * dr.x, cf_bond * dr.y, cf_bond * dr.z);
  
  d_virial.x += dr.x * Fbond.x;
  d_virial.y += dr.y * Fbond.y;
  d_virial.z += dr.z * Fbond.z;

  lap_conf   += itrs.cf_spring * (6.0 * Parameter::i_bleng - 4.0 * inv_dr);
  
  force[l_dst[0]] -= Fbond;
  force[l_dst[1]] += Fbond;
}

template<>
void F_calculator::CalcBondBend<1>(const double3* __restrict	pr,
				   double3*       __restrict    force,
				   double3&       __restrict    d_virial,
				   double&        __restrict    lap_conf,			   
				   const int			beg_idx,
				   const int*			elem_idx,
				   const Parameter&		param)
{}

F_calculator::F_calculator(const Parameter& param) {
  itrs  = param.GetIntractions();
  binfo = param.GetBindInform();
  const int th_numb   = omp_get_max_threads();
  buf_vir.resize(th_numb, double3(0.0, 0.0, 0.0));
  buf_lap_pot.resize(th_numb, 0.0);
  DevideCell(param);
}

void F_calculator::AddConservForce(const double3* __restrict	pr,
				   const par_prop* __restrict   prop,
				   double3* __restrict		force,
				   const CellList&		clist,
				   const Parameter&	        param,
				   ChemInfo&			cheminfo) {
  ClearForce(force);

#pragma omp parallel 
  {
    CalcForceHalf(pr, prop, force, clist, param, 0);
#pragma omp barrier
    CalcForceHalf(pr, prop, force, clist, param, 1);
  }
  
  CalcForceInAmp(pr,force,cheminfo,param);
}

void F_calculator::AddDisspRandom(const double3*  __restrict pr,
				  double3*        __restrict pv,
				  const par_prop* __restrict prop,
				  const CellList&            clist,
				  const Parameter&           param,
				  RNG&                       rng) {
#pragma omp parallel
  {
    AddCollisionHalf(pr, pv, prop, clist, param, rng, 0);
#pragma omp barrier
    AddCollisionHalf(pr, pv, prop, clist, param, rng, 1);
  }
}

//ASSUME: tailN <= headN
void F_calculator::RegistNearIdx(const int&	itr_id,
				 const int&     tar_id,
				 const float&	dist,
				 ChemInfo&	cheminfo) {
  const bool chem_flag = (cheminfo.prtcl_chem[itr_id] == false) && (cheminfo.prtcl_chem[tar_id] == false);
  if (chem_flag) {
    const int i_unit = cheminfo.lipid_unit[itr_id];
    const int t_unit = cheminfo.lipid_unit[tar_id];
      
    const bool t_tail_i_head = (t_unit == Parameter::REAC_PART) && (i_unit == Parameter::REAC_PART-1);
    const bool t_head_i_tail = (t_unit == Parameter::REAC_PART-1) && (i_unit == Parameter::REAC_PART);
    if (t_head_i_tail) {
      if (dist < cheminfo.near_info[itr_id].dist) {
	cheminfo.near_info[itr_id].dist = dist;
	cheminfo.near_info[itr_id].idx  = tar_id;
      }
      return;
    }
    if (t_tail_i_head) {
      if (dist < cheminfo.near_info[tar_id].dist) {
	cheminfo.near_info[tar_id].dist = dist;
	cheminfo.near_info[tar_id].idx  = itr_id;
      }
      return;
    }
  }
}

void F_calculator::CountNearWater(const int& itr_id,
				  const int& tar_id,		      
				  const par_prop& i_prp,
				  const par_prop& t_prp,
				  float dist,
				  ChemInfo& cheminfo) {
  if (cheminfo.lipid_unit[itr_id] == Parameter::REAC_PART) {
    const bool t_is_water = (t_prp == Water);  
    const bool i_chm_flag = cheminfo.prtcl_chem[itr_id];      
    if (t_is_water && i_chm_flag) cheminfo.near_water[itr_id] += (dist < Parameter::ch_leng);
    //cheminfo.near_water[itr_id] += 1.0 / (expf(20.0*(dist-0.98)) + 1.0);
    return;
  } else if (cheminfo.lipid_unit[tar_id] == Parameter::REAC_PART) {
    const bool i_is_water = (i_prp == Water);      
    const bool t_chm_flag = cheminfo.prtcl_chem[tar_id];
    if (i_is_water && t_chm_flag) cheminfo.near_water[tar_id] += (dist < Parameter::ch_leng);
      //cheminfo.near_water[tar_id] += 1.0 / (expf(20.0*(dist-0.98)) + 1.0);
    return;
  }
}

void F_calculator::CalcForceHalf(const double3*  __restrict pr,
				 const par_prop* __restrict prop,
				 double3*        __restrict force,
				 const CellList&            clist,
				 const Parameter&           param,
				 const int                  call_num) {
  double3 sum_vir(0.0);     double d_lap_pot(0.0);
  const int	tid	 = omp_get_thread_num();
  const int	beg_grid = grid_numbtw * (tid * numb_band + p_num_band[call_num]) ;
  const int	end_grid = beg_grid + grid_numbtw * p_num_band[call_num + 1];
  const int	begi	 = clist.buck_addrs[beg_grid];
  const int	endi	 = clist.buck_addrs[end_grid];
  int		begj	 = 0;
  int		bef_grid = -1;
  for (int i = begi; i < endi; i++) {
    double3		sumF_a(0.0);
    const int		pi	 = clist.prtcl_idx[i];
    const int		tar_grid = clist.prtcl_in_cell[pi];
    const double3	ri       = pr[pi];
    const par_prop	prp_pi	 = prop[pi];
    if (tar_grid != bef_grid) {
      begj = 0;
      bef_grid = tar_grid;
    }
    begj++;
    
    const int endj = clist.n_near_prtcl[tar_grid];
    for (int j = begj; j < endj; j++) {
      const int pj = clist.near_prtcl_idx[tar_grid][j];
      double3 dr(ri.x - pr[pj].x, ri.y - pr[pj].y, ri.z - pr[pj].z);
      MinImage(dr, param);
      const double dr2  = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
      if (dr2 < 1.0) {
#ifdef DEBUG
	/*#pragma omp critical
	{
	  if(pi < pj){std::cout << pi << " " << pj << std::endl;
	  }else{std::cout << pj << " " << pi << std::endl;}
	  }*/
#endif
	const double inv_dr   = 1.0 / std::sqrt(dr2);
	const par_prop prp_pj = prop[pj];
	const double dr_m_one = inv_dr - 1.0;
	const double  cf_c    = itrs.cf_repul[prp_pi][prp_pj];
	const double all_cf   = cf_c * dr_m_one;
	const double3 dF_a(all_cf * dr.x, all_cf * dr.y, all_cf * dr.z);

	sumF_a    += dF_a;
	force[pj] -= dF_a;
	
	sum_vir.x += dF_a.x*dr.x;
	sum_vir.y += dF_a.y*dr.y;
	sum_vir.z += dF_a.z*dr.z;
	
	d_lap_pot += cf_c * (6.0 - 4.0 * inv_dr);
	
#ifdef CHEM_MODE	
	const float dist = dr2*inv_dr;
	RegistNearIdx(pi,pj,dist,cheminfo);
#ifdef INVERSE_CHEM
	CountNearWater(pi,pj,prp_pi,prp_pj,dist,cheminfo);
#endif
#endif
	
      } //end if
    } //end j loop
    force[pi] += sumF_a;
  } //end i loop  
  
  buf_vir[tid]     += sum_vir;
  buf_lap_pot[tid] += d_lap_pot;
}

void F_calculator::AddCollisionHalf(const double3*  __restrict pr,
				    double3*        __restrict pv,
				    const par_prop* __restrict prop,
				    const CellList&            clist,
				    const Parameter&           param,
				    RNG&                       rng,
				    const int                  call_num) {
  const double nois_amp  = std::sqrt(2.0 * param.tempera);
  
  const int	tid	 = omp_get_thread_num();
  const int	beg_grid = grid_numbtw * (tid * numb_band + p_num_band[call_num]) ;
  const int	end_grid = beg_grid + grid_numbtw * p_num_band[call_num + 1];
  const int	begi	 = clist.buck_addrs[beg_grid];
  const int	endi	 = clist.buck_addrs[end_grid];
  int		begj	 = 0;
  int		bef_grid = -1;
  for (int i = begi; i < endi; i++) {
    const int		pi	 = clist.prtcl_idx[i];
    const int		tar_grid = clist.prtcl_in_cell[pi];
    const double3	ri	 = pr[pi];
    const par_prop	prp_pi	 = prop[pi];
    if (tar_grid != bef_grid) {
      begj=0;
      bef_grid=tar_grid;
    }
    begj++;
    const int endj = clist.n_near_prtcl[tar_grid];
    
    for (int j = begj; j < endj; j++) {
      const int pj = clist.near_prtcl_idx[tar_grid][j];
      double3 dr(ri.x - pr[pj].x, ri.y - pr[pj].y, ri.z - pr[pj].z);
      MinImage(dr, param);
      const double dr2  = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
      if (dr2 < 1.0) {
	const double    dr_norm = std::sqrt(dr2);
	const par_prop	prp_pj  = prop[pj];
	const double	nrml    = rng.Normal(tid) * nois_amp;
	const double	inv_dr  = 1.0 / dr_norm;

#ifdef VISC_KERN_SQRT
	Collision_pow_half(pv[pi], pv[pj], dr, itrs.cf_gamma[prp_pi][prp_pj], nrml, inv_dr, dr_norm, param);
#else
	Collision_pow_one(pv[pi], pv[pj], dr, itrs.cf_gamma[prp_pi][prp_pj], nrml, inv_dr, param);
#endif

      } //end if
    } //end j loop
  } //end i loop
}

void F_calculator::CalcForceInAmp(const double3* __restrict	pr,
				  double3* __restrict		force,
				  const ChemInfo&		cheminfo,
				  const Parameter&		param) {
  double3 d_virial(0.0); double d_lap_potent = 0.0;
  for (int i = 0; i < Parameter::sys_size; i++) {
    if (cheminfo.lipid_unit[i] == Parameter::REAC_PART) {
      if (cheminfo.prtcl_chem[i]) {
	const int beg_idx = cheminfo.lipid_idx[i] * Parameter::ALL_UNIT_N;
	CalcBondBend<Parameter::ALL_UNIT_N>(pr, force,d_virial, d_lap_potent, beg_idx, cheminfo.lip_elem_idx, param);
      } else {
	const int beg_idx = cheminfo.part_idx[i] * Parameter::TAIL_PART;
	CalcBondBend<Parameter::TAIL_PART >(pr, force,d_virial, d_lap_potent, beg_idx, cheminfo.tail_elem_idx,param);
      }
    } else if (cheminfo.lipid_unit[i] == Parameter::REAC_PART - 1) {
      if (!cheminfo.prtcl_chem[i]) {
	const int beg_idx = cheminfo.part_idx[i] * Parameter::REAC_PART;
	CalcBondBend<Parameter::REAC_PART>(pr,force,d_virial,d_lap_potent,beg_idx,cheminfo.head_elem_idx,param);
      }
    }
  }
  
  buf_vir[0] += d_virial;
  buf_lap_pot[0] += d_lap_potent;
}

//NOTE: membrane is placed parallel to xz plane.
void F_calculator::AddBindForceCyl(const double3*  __restrict pr,
				   const par_prop* __restrict prop,
				   double3*        __restrict force,
				   const ChemInfo&  cheminfo) {
  const double bind_cutof = 2.0;
  for (int i = 0; i < Parameter::sys_size; i++) {
    const int l_idx = cheminfo.lipid_idx[i];
    bool chem = false;
    if (l_idx != -1) chem = cheminfo.lipid_chem[l_idx];
    if ((!chem) && (prop[i] == Hyphob)) {
      const double3 dr(pr[i].x - binfo.bind_center.x,
		       0.0,
		       pr[i].z - binfo.bind_center.z);
      const double dist = dr.norm2();
      const double bound2prtcl = dist - binfo.bind_radius;
      if ((bound2prtcl > 0.0) && (bound2prtcl < bind_cutof)) {
	const double inv_dist = 1.0 / dist;
	const double cf_bind = binfo.bind_coef*bound2prtcl*inv_dist;
	const double3 bF(cf_bind*dr.x,
			 0.0,
			 cf_bind*dr.z);
	force[i] -= bF;
      }
    }
  }
}

void F_calculator::AddBindForceSph(const double3* __restrict pr,
				   const par_prop* __restrict prop,
				   double3* __restrict force) {
  const double bind_cutoff = Parameter::ALL_UNIT_N * Parameter::b_leng * 0.5;
  for (int i = 0; i < Parameter::sys_size; i++) {
    if (prop[i] != Water) {
      const double3 dr(pr[i].x - binfo.bind_center.x,
		       pr[i].y - binfo.bind_center.y,
		       pr[i].z - binfo.bind_center.z);
      const double dist = dr.norm2();
      if (fabs(binfo.bind_radius - dist) < bind_cutoff) {
	const double dr_m_one = binfo.bind_radius / dist - 1.0;
	const double cf_bind = binfo.bind_coef * dr_m_one;
	const double3 bF(cf_bind * dr.x,
			 cf_bind * dr.y,
			 cf_bind * dr.z);
	force[i] += bF;
      }
    }
  }
}

void F_calculator::AddPosRes(const double3* __restrict pr,
			     const double3* __restrict pr_base,
			     double3* __restrict force,
			     const bool* __restrict rest_on) {
  constexpr double bcf = 200.0;
  
  for (int i = 0; i < Parameter::sys_size; i++) {
    if (rest_on[i]) {
      const double3 dr(pr[i].x - pr_base[i].x,
		       0.0, //pr[i].y - pr_base[i].y,
		       pr[i].z - pr_base[i].z);
      const double3 bF(bcf * dr.x,
		       0.0, //bcf * dr.y,
		       bcf * dr.z);
      force[i] -= bF;
    }
  }
}
