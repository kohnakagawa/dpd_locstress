#pragma once

#include "parameter.hpp"
#include "rng.hpp"

#define INLINE __attribute__((always_inline))

class ChemManager{
  INLINE bool CalcChemProb(const float& dist,
			   const Parameter& param,
			   RNG& rng)
  {
    const float t_prob = rng.Uniform(0);
    const float c_prob = (param.coef_prob[0]) * (dist < param.coef_prob[1]);
    const bool ret = (t_prob < c_prob);
    return ret;
  }

  void DetChemReac(const Parameter& param, ChemInfo& cheminfo, RNG& rng);
  void ReSearchNearest(ChemInfo& cheminfo);
  
public:
  INLINE void ChemEvent(const Parameter& param, ChemInfo& cheminfo, RNG& rng){
    ReSearchNearest(cheminfo);
    DetChemReac(param,cheminfo,rng);
  }
  void ClearNearInfo(ChemInfo& cheminfo);
  
  //NOTE:this function should be called after bucket sort
  void RegistLipidIdx(ChemInfo& cheminfo);
};

#undef INLINE
