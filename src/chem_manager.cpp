#include "chem_manager.hpp"

void ChemManager:: DetChemReac(const Parameter& param,
			       ChemInfo& cheminfo,
			       RNG& rng)
{
  //ASSUME: tailN <= headN
  for(int i=0; i<Parameter::SYS_SIZE; i++){
    const int   unit_id = cheminfo.lipid_unit[i];
    const bool chm_flag = cheminfo.prtcl_chem[i];
#ifdef INVERSE_CHEM
    if(unit_id == param.REAC_PART){
      if(chm_flag){
	//divided
	const int l_idx    = cheminfo.lipid_idx[i];
	//const float weight = (cheminfo.near_water[i] - param.prob_cutof) * 0.077;
	const float weight = (cheminfo.near_water[i] != 0.0);
	const float t_prob = rng.Uniform(0);
	const float c_prob = param.coef_prob[2] * weight;
	if(t_prob < c_prob){
	  cheminfo.lipid_chem[l_idx] = false;
	}
      }
    }else
#endif
      if(unit_id == param.REAC_PART-1){
	if(!chm_flag){
	  //connected
	  const float_int n_info = cheminfo.near_info[i];
	  const bool   chem_crit = CalcChemProb(n_info.dist,param,rng);
	  if(chem_crit && (n_info.idx != -1)){
	    cheminfo.prtcl_chem[i         ] = true;
	    cheminfo.prtcl_chem[n_info.idx] = true;
	    const int l_idx = cheminfo.lipid_idx[n_info.idx];
	    cheminfo.lipid_idx[i] = l_idx;
	    cheminfo.lipid_chem[l_idx]     = true;
	    
	    //register lip_elem_idx
	    const int base_id = Parameter::ALL_UNIT_N * l_idx;
	    const int part_id[2] = {cheminfo.part_idx[i],cheminfo.part_idx[n_info.idx]};
	  
	    for(int unit=0; unit<Parameter::REAC_PART; unit++){
	      const int h_e_idx = cheminfo.head_elem_idx[part_id[0]*Parameter::REAC_PART + unit];
	      cheminfo.lip_elem_idx[base_id + unit] = h_e_idx;
	      cheminfo.lipid_idx[h_e_idx] = l_idx;
	    }
	    for(int unit=Parameter::REAC_PART; unit<Parameter::ALL_UNIT_N; unit++){
	      const int t_e_idx = cheminfo.tail_elem_idx[part_id[1]*Parameter::TAIL_PART + unit - Parameter::REAC_PART];
	      cheminfo.lip_elem_idx[base_id + unit] = t_e_idx;
	      cheminfo.lipid_idx[t_e_idx] = l_idx;
	    }
	  }
	}
      }
  }

  //ASSUME: tailN <= headN
  for(int i=0; i<Parameter::SYS_SIZE; i++){
    bool flag = true;
    const int l_idx = cheminfo.lipid_idx[i    ];
    if(l_idx != -1) flag = cheminfo.lipid_chem[l_idx];
    if(!flag){
      cheminfo.prtcl_chem[i] = false;
      const int unit = cheminfo.lipid_unit[i];
      if(unit < param.REAC_PART) cheminfo.lipid_idx[i] = -1;
    }
  }
}

//ASSUME: tailN <= headN
void ChemManager::ReSearchNearest(ChemInfo& cheminfo) {
  for(int tail_idx=0; tail_idx<Parameter::SYS_SIZE; tail_idx++){
    const int head_idx = cheminfo.near_info[tail_idx].idx;
    if( (head_idx != -1) && (cheminfo.lipid_unit[tail_idx] == Parameter::REAC_PART)){
      const float tail2head = cheminfo.near_info[tail_idx].dist;
      cheminfo.near_info[tail_idx].dist = 4.0;
      cheminfo.near_info[tail_idx].idx  = -1;
      if(cheminfo.near_info[head_idx].dist > tail2head){
	cheminfo.near_info[head_idx].dist = tail2head;
	cheminfo.near_info[head_idx].idx  = tail_idx;
      }
    }
  }
}

void ChemManager::ClearNearInfo(ChemInfo& cheminfo)
{
  for(int i=0; i<Parameter::SYS_SIZE; i++){
    cheminfo.near_info[i].dist = 4.0;
    cheminfo.near_info[i].idx  = -1;
    cheminfo.near_water[i] = 0.0;
  }
}

void ChemManager::RegistLipidIdx(ChemInfo& cheminfo) {
  for(int i=0; i<Parameter::SYS_SIZE; i++){
    bool chm_flag = false;
    const int l_idx     = cheminfo.lipid_idx[i];
    const int unit_id   = cheminfo.lipid_unit[i];
    if(l_idx != -1) chm_flag = cheminfo.lipid_chem[l_idx];
    
    if(chm_flag){
      cheminfo.lip_elem_idx[Parameter::ALL_UNIT_N*l_idx + unit_id] = i;
    }
    
    const int p_idx = cheminfo.part_idx[i];
    if(p_idx != -1){
      if(unit_id < Parameter::REAC_PART){
	cheminfo.head_elem_idx[Parameter::REAC_PART*p_idx + unit_id] = i;      
      }else{
	cheminfo.tail_elem_idx[Parameter::TAIL_PART*p_idx + unit_id - Parameter::REAC_PART] = i;
      }
    }
  }
}
