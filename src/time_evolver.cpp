#include "dpdsystem.hpp"

void dpdsystem::VelocVerlet1() {
#pragma omp parallel for
  for (int i = 0; i < Parameter::sys_size; i++) {
    force_bef[i] = force[i];
    pv_bef[i]    = pv[i];
    const double3 temp_dr((pv[i].x + 0.5 * p_param->dt * force[i].x) * p_param->dt,
			  (pv[i].y + 0.5 * p_param->dt * force[i].y) * p_param->dt,
			  (pv[i].z + 0.5 * p_param->dt * force[i].z) * p_param->dt);
    
    delta_sumr[i] += temp_dr;
    
    const double3 temp_r(pr[i].x + temp_dr.x, 
			 pr[i].y + temp_dr.y, 
			 pr[i].z + temp_dr.z);

    pr[i].x = temp_r.x - std::floor(temp_r.x * p_param->iL.x) * p_param->L.x;

#ifdef Y_REFLECT_BOUND
    ReflectBound(temp_r.y, pr[i].y, pv_bef[i], force_bef[i], p_param->L.y);
#else
    pr[i].y = temp_r.y - std::floor(temp_r.y * p_param->iL.y) * p_param->L.y;
#endif    
    
    pr[i].z = temp_r.z - std::floor(temp_r.z * p_param->iL.z) * p_param->L.z;
    
    pv[i].x = pv_bef[i].x + 0.5 * p_param->dt * force_bef[i].x;
    pv[i].y = pv_bef[i].y + 0.5 * p_param->dt * force_bef[i].y;
    pv[i].z = pv_bef[i].z + 0.5 * p_param->dt * force_bef[i].z;
  }
}

void dpdsystem::VelocVerlet2(){
  for (int i = 0; i < Parameter::sys_size; i++) {
    pv[i].x = pv_bef[i].x + 0.5 * (force[i].x + force_bef[i].x) * p_param->dt;
    pv[i].y = pv_bef[i].y + 0.5 * (force[i].y + force_bef[i].y) * p_param->dt;
    pv[i].z = pv_bef[i].z + 0.5 * (force[i].z + force_bef[i].z) * p_param->dt;
  }
}
