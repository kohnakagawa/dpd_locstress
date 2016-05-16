#pragma once

#include "dpdsystem.hpp"
#include "rng.hpp"
#include <utility>

class Initializer
{
  int h_array_id, t_array_id, amp_array_id;

  void InitVeloc(dpdsystem &sDPD,const Parameter& param,RNG& rng);
  void ReviseVeloc(dpdsystem &sDPD,const Parameter& param)const;

  void InitConfig(dpdsystem &sDPD,const Parameter& param,RNG& rng);
  void ApplyPeriodicBound(double3& r,const Parameter& param);
  void SetPartialLipidId(dpdsystem& sDPD, const int i, int& base_idx, const int unit_leng, const int elem);
  
  //NOTE: r:base position n:membrane's normal vector
  //head case offset=0                 end=param.head_part_n
  //tail case offset=param.head_part_n end=param.tail_part_n
  template<int offset,int end>
  void SetAmphilPartPos(const double3& r, const double3& n,const Parameter& param,int& phil_idx,int& oil_idx, dpdsystem& sDPD,bool is_head){
    for(int unit=0; unit<end; unit++){
      double3 temp_pos = r + (unit + offset) * param.b_leng * n;
      ApplyPeriodicBound(temp_pos,param);

        //hydrophilic part
      if(unit + offset < Parameter::HYPHIL_N){
	sDPD.SetPosition(phil_idx,temp_pos);
	sDPD.SetProp(phil_idx,Hyphil);
	sDPD.SetLipidUnit(phil_idx,unit+offset);
	
	if(is_head){
	  sDPD.SetHeadElemIdx(h_array_id,phil_idx);
	  sDPD.SetPartIdx(phil_idx,h_array_id / Parameter::REAC_PART);
	  h_array_id++;
	}else{
	  sDPD.SetTailElemIdx(t_array_id,phil_idx);
	  sDPD.SetPartIdx(phil_idx,t_array_id / Parameter::TAIL_PART);
	  t_array_id++;
	}
	if(amp_array_id < param.ampN * Parameter::ALL_UNIT_N){
	  sDPD.SetLipidElemIdx(amp_array_id,phil_idx);
	  amp_array_id++;
	}
	
	phil_idx++;
	
	//hydrophobic part
      }else{
	sDPD.SetPosition(oil_idx,temp_pos);
	sDPD.SetProp(oil_idx,Hyphob);
	sDPD.SetLipidUnit(oil_idx,unit+offset);
	
	if(is_head){
	  sDPD.SetHeadElemIdx(h_array_id,oil_idx);
	  sDPD.SetPartIdx(oil_idx,h_array_id / Parameter::REAC_PART);
	  h_array_id++;
	}else{
	  sDPD.SetTailElemIdx(t_array_id,oil_idx);
	  sDPD.SetPartIdx(oil_idx,t_array_id / Parameter::TAIL_PART);
	  t_array_id++;
	}
	if(amp_array_id < Parameter::ALL_UNIT_N*param.ampN){
	  sDPD.SetLipidElemIdx(amp_array_id,oil_idx);
	  amp_array_id++;
	}
	
	oil_idx++;
      }
    }
  }
  
  void ClearLipidIdx(dpdsystem& sDPD);
  void SetLipidId(dpdsystem &sDPD, const Parameter& param, int elem);
  void MakeRandomDispert(int elem,int& phil_idx,int& oil_idx,const Parameter& param,dpdsystem &sDPD,RNG& rng);
  void MakeDisk(int elem, int& phil_idx, int& oil_idx, const double3& offset,
		const Parameter &param, dpdsystem &sDPD);
  void MakeSheet(int elem,int& phil_idx,int& oil_idx,double x_leng, double y_leng, double z_leng,
		 const Parameter &param,dpdsystem &sDPD,RNG& rng);
  void MakeHalfSheet(int elem,int& phil_idx,int& oil_idx,double x_leng, double y_leng, double z_leng, const Parameter &param,dpdsystem &sDPD,RNG& rng,bool flag);
  double MakeCylindSheet(int& phil_idx, int& oil_idx, double cylrad, double cyl_Ly, double offset, int cylelem, const Parameter &param, dpdsystem &sDPD, RNG& rng);
  
  void MakeArcLine(int elem, double d_the, double3 center, double rad, double sign, int& phil_idx, int& oil_idx, dpdsystem &sDPD, const Parameter &param);
  
  bool MakeSphLine(const double the, const double3& center, const double rad, const double sign, const double lip_len,
		   const int end_phil_idx, int& phil_idx, int& oil_idx, dpdsystem& sDPD, const Parameter& param);

  void MakeForEachTheta(const int elem, const double d_the, const double offset, const double rad, const double sign, 
			const double lip_len, const int end_phil_idx, int& phil_idx, int& oil_idx, dpdsystem& sDPD, 
			const Parameter& param);

  void MakeSphereSheet(int elem, int& phil_idx, int&oil_idx, const Parameter& param, dpdsystem& sDPD);
  
  //NOTE: These functions are applicable for rectangular bilayer sheet.
  void MakeWater(int elem, int& w_idx, const double3& reg_dw, const double3& reg_up, const Parameter &param, dpdsystem &sDPD, RNG& rng, bool flag);
  void MakeWaterAsym(int elem, int& w_idx, const double3& reg_dw, const double3& reg_up, const Parameter &param, dpdsystem &sDPD, RNG& rng, const double ratio);
  void MakeWaterCyl(int elem, int& w_idx, double cyl_rad, double3 cent, const double* y_reg, const Parameter &param, dpdsystem &sDPD, RNG& rng);
  void MakeWaterSphere(const int elem, int& w_idx, int in_elem, const Parameter& param, dpdsystem& sDPD, RNG& rng);
  void MakeEmbedTail(int elem, int& phil_idx, int& oil_idx, const double3& reg_dw, const double3& reg_up, const Parameter &param,dpdsystem &sDPD,RNG& rng);
  void MakeHeadInWater(int elem, int& phil_idx, int& oil_idx, const double3& reg_dw, const double3& reg_up, const Parameter &param, dpdsystem &sDPD,RNG& rng);
  void MakeHeadInWaterCyl(int elem, int& phil_idx, int& oil_idx, double cyl_rad, double3 cent, double z_hei, bool flag, const Parameter &param, dpdsystem &sDPD, RNG& rng);
  void MakeHeadInWaterSph(int elem, int& phil_idx, int& oil_idx, const Parameter &param, dpdsystem &sDPD, RNG& rng);
  void MakeOilbulk(int elem, int& phil_idx, int& oil_idx, const double3& reg_dw, const double3& reg_up, const Parameter& param, dpdsystem& sDPD, RNG& rng);
  
  enum{
    NORMAL = 0,      //hN == bN
    TAIL_EMBEDDED,  //hN <  bN
    HEAD_DISPERT,//hN >  bN
  };
  
  void GenDiskLam(dpdsystem &sDPD,const Parameter& param,RNG& rng,int mode);
  void GenPlaneLam(dpdsystem &sDPD,const Parameter& param,RNG& rng,int mode);
  void GenAsynPlaneLam(dpdsystem &sDPD,const Parameter& param,RNG& rng,int mode);
  void GenStrip(dpdsystem &sDPD,const Parameter& param, RNG& rng,int mode);
  void GenRandomSolute(dpdsystem &sDPD, const Parameter& param,RNG& rng);
  void GenOilDrop(dpdsystem &sDPD, const Parameter& param, RNG& rng);
  void GenOilCuboid(dpdsystem &sDPD, const Parameter& param, RNG& rng);
  void GenCylindLam(dpdsystem &sDPD, const Parameter& param, RNG& rng, int mode);
  void GenSphere(dpdsystem &sDPD, const Parameter& param, RNG& rng);
  
  void LoadRestartConfig(dpdsystem& sDPD, ChemInfo& cheminfo, const Parameter& param, std::ifstream& fin) const;
  void NumberIsMatched(const int r_v, const int l_v, const std::string& mess) const ;
  void SetBondedParameter(dpdsystem& sDPD, ChemInfo& cheminfo, const Parameter& param);
  void CheckRestartConfigIsValid(const dpdsystem& sDPD, const ChemInfo& cheminfo) const;

  void CheckInited(const dpdsystem &sDPD, const Parameter& param) const;
  
public:
  Initializer():h_array_id(0), t_array_id(0), amp_array_id(0){};
  void GenParticles(dpdsystem &sDPD, ChemManager& chemmanage, ChemInfo& cheminfo,const Parameter& param, RNG& rng);
};
