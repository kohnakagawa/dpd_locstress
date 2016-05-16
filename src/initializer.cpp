#include "initializer.hpp"
#include "chem_manager.hpp"

void Initializer::InitVeloc(dpdsystem &sDPD,const Parameter& param,RNG& rng){
  for(int i=0; i<Parameter::SYS_SIZE; i++){
    const double3 tempv(rng.Normal(0,0.0,sqrt(param.tempera)),
			rng.Normal(0,0.0,sqrt(param.tempera)),
			rng.Normal(0,0.0,sqrt(param.tempera)));
    sDPD.SetMoment(i,tempv);
  }
  ReviseVeloc(sDPD,param);
}

void Initializer::ReviseVeloc(dpdsystem &sDPD,const Parameter& param)const{
  double3 vsum(0.0);
  const double  invN = 1.0/Parameter::SYS_SIZE;

  const double3* v = sDPD.pv;
  for(int i=0; i<Parameter::SYS_SIZE; i++){
    vsum += v[i];
  }

  vsum *= invN;

  for(int i=0; i<Parameter::SYS_SIZE; i++){
    double3 tempv = v[i];
    tempv -= vsum;
    sDPD.SetMoment(i,tempv);
  }
}

void Initializer::InitConfig(dpdsystem &sDPD,const Parameter& param,RNG& rng){
  //GenDiskLam(sDPD,param,rng,NORMAL);
  //GenDiskLam(sDPD,param,rng,TAIL_EMBEDDED);
  //GenDiskLam(sDPD,param,rng,HEAD_DISPERT);

  //GenPlaneLam(sDPD,param,rng,NORMAL);
  //GenPlaneLam(sDPD,param,rng,TAIL_EMBEDDED);
  //GenPlaneLam(sDPD,param,rng,HEAD_DISPERT);
  //GenAsynPlaneLam(sDPD,param,rng,NORMAL);
  GenAsynPlaneLam(sDPD,param,rng,TAIL_EMBEDDED);
    
  //GenStrip(sDPD,param,rng,NORMAL);
  //GenStrip(sDPD,param,rng,TAIL_EMBEDDED);
    
  //GenCylindLam(sDPD,param,rng,NORMAL);
  //GenCylindLam(sDPD,param,rng,HEAD_DISPERT);

  //GenSphere(sDPD, param, rng);

  //GenRandomSolute(sDPD,param,rng);

  //GenOilDrop(sDPD,param,rng);

  //GenOilCuboid(sDPD,param,rng);
}

void Initializer::ApplyPeriodicBound(double3& r,const Parameter& param){
  r.x -= floor(r.x * param.iL.x) * param.L.x;
  r.y -= floor(r.y * param.iL.y) * param.L.y;
  r.z -= floor(r.z * param.iL.z) * param.L.z;
}

void Initializer::SetPartialLipidId(dpdsystem& sDPD, const int i, int& base_idx, const int unit_leng, const int elem){
  for(int unit=0; unit<unit_leng; unit++){
    sDPD.SetLipidIdx(base_idx, i);
    if(i < elem) sDPD.SetChemConf(base_idx, true);
    base_idx++;
  }
}

void Initializer::ClearLipidIdx(dpdsystem& sDPD){
  for(int i=0; i<Parameter::SYS_SIZE; i++){
    sDPD.SetLipidIdx(i, -1);
    sDPD.SetPartIdx(i, -1);
    sDPD.SetLipidUnit(i, -1);
    sDPD.SetChemConf(i, false);
  }
}

void Initializer::SetLipidId(dpdsystem &sDPD,const Parameter& param, int elem){
  assert(elem <= param.ampN);
  
  //hydrophilic: base_idx[0] hydrophobic: base_idx[1]
  int base_idx[2] = {0, param.hN};
  
  for(int i=0; i<elem; i++){
    SetPartialLipidId(sDPD, i, base_idx[0], Parameter::HYPHIL_N, elem);
    SetPartialLipidId(sDPD, i, base_idx[1], Parameter::HYPHOB_N, elem);
    sDPD.SetLipidChem(i, true);
  }

  if(Parameter::HYPHIL_N >= Parameter::REAC_PART){
    if(param.tailN <= param.headN) base_idx[0] += Parameter::REAC_PART * (param.headN - elem);
  }else{
    if(param.tailN <= param.headN) base_idx[1] += (Parameter::REAC_PART - Parameter::HYPHIL_N) * (param.headN - elem);
  }
  
  for(int i=elem; i<param.ampN; i++){
    if(Parameter::HYPHIL_N >= Parameter::REAC_PART){
      if(param.tailN <= param.headN){
	//tail
	SetPartialLipidId(sDPD, i, base_idx[0], Parameter::HYPHIL_N - Parameter::REAC_PART, elem);
	SetPartialLipidId(sDPD, i, base_idx[1], Parameter::HYPHOB_N, elem);
      }else{
	//head
	SetPartialLipidId(sDPD, i, base_idx[0], Parameter::REAC_PART, elem);
      }
    }else{
      if(param.tailN <= param.headN){
	//tail
	SetPartialLipidId(sDPD, i, base_idx[1], Parameter::TAIL_PART, elem);
      }else{
	//head
	SetPartialLipidId(sDPD, i, base_idx[0], Parameter::HYPHIL_N, elem);
	SetPartialLipidId(sDPD, i, base_idx[1], Parameter::REAC_PART - Parameter::HYPHIL_N, elem);
      }
    }
    sDPD.SetLipidChem(i, false);
  }
}

void Initializer::MakeRandomDispert(int elem,int& phil_idx,int& oil_idx,const Parameter& param,dpdsystem &sDPD,RNG& rng){
  const int philic_end_idx = phil_idx + elem;
  while(phil_idx < philic_end_idx){
    double3 n_v(rng.Uniform(0),
		rng.Uniform(0),
		rng.Uniform(0));
    n_v /= n_v.norm2();
      
    const double3 r(rng.Uniform(0) * param.L.x,
		    rng.Uniform(0) * param.L.y,
		    rng.Uniform(0) * param.L.z);
    
    SetAmphilPartPos<0,Parameter::REAC_PART>(r,n_v,param,phil_idx,oil_idx,sDPD,true);
    SetAmphilPartPos<Parameter::REAC_PART,Parameter::TAIL_PART>(r,n_v,param,phil_idx,oil_idx,sDPD,false);
  }
}

void Initializer::MakeDisk(int elem, int& phil_idx, int& oil_idx, const double3& offset,
			   const Parameter &param,dpdsystem &sDPD)
{
  std::pair<int,int> ii(0,0);
  const double lat_leng = sqrt(param.lip_area);
  const int philic_end_idx = phil_idx + elem;
  
  const double up_x = param.L.x - offset.x;
  
  while(phil_idx < philic_end_idx){
    const double3 up(ii.first * lat_leng + offset.x, param.hL.y + Parameter::ALL_UNIT_N*param.b_leng + 0.01 + offset.y , ii.second * lat_leng + offset.z);
    const double3 dw(ii.first * lat_leng + offset.x, param.hL.y - Parameter::ALL_UNIT_N*param.b_leng - 0.01 + offset.y , ii.second * lat_leng + offset.z);
    
    double3 n_v(0.0,1.0,0.0);
    
    SetAmphilPartPos<0,Parameter::REAC_PART>(dw,n_v,param,phil_idx,oil_idx,sDPD,true);
    SetAmphilPartPos<Parameter::REAC_PART,Parameter::TAIL_PART>(dw,n_v,param,phil_idx,oil_idx,sDPD,false);
    
    if(phil_idx >= philic_end_idx){
      break;
    }
    
    n_v.y *= -1.0;

    SetAmphilPartPos<0,Parameter::REAC_PART>(up,n_v,param,phil_idx,oil_idx,sDPD,true);
    SetAmphilPartPos<Parameter::REAC_PART,Parameter::TAIL_PART>(up,n_v,param,phil_idx,oil_idx,sDPD,false);
      
    ii.first++;
    if(ii.first * lat_leng + offset.x > up_x){
      ii.first = 0;
      ii.second++;
    }
    if(ii.second * lat_leng + offset.z >= param.L.z) std::cout << ii.second * lat_leng + offset.z << std::endl;
    assert(ii.second * lat_leng + offset.z < param.L.z);
  }
}

void Initializer::MakeSheet(int elem,int& phil_idx,int& oil_idx,double x_leng, double y_leng, double z_leng,
			    const Parameter &param,dpdsystem &sDPD,RNG& rng)
{
  const int philic_end_idx = phil_idx + elem;
    
  while(phil_idx < philic_end_idx){
    const double3 up(x_leng * rng.Uniform(0),
		     y_leng + Parameter::ALL_UNIT_N*param.b_leng + 0.01,
		     z_leng * rng.Uniform(0));
    const double3 dw(x_leng * rng.Uniform(0),
		     y_leng - Parameter::ALL_UNIT_N*param.b_leng - 0.01,
		     z_leng * rng.Uniform(0));
    double3 n_v(0.0,1.0,0.0);
    
    SetAmphilPartPos<0,Parameter::REAC_PART>(dw,n_v,param,phil_idx,oil_idx,sDPD,true);
    SetAmphilPartPos<Parameter::REAC_PART,Parameter::TAIL_PART>(dw,n_v,param,phil_idx,oil_idx,sDPD,false);
    
    n_v.y *= -1.0;
    if(phil_idx >= philic_end_idx){
      break;
    }
    
    SetAmphilPartPos<0,Parameter::REAC_PART>(up,n_v,param,phil_idx,oil_idx,sDPD,true);
    SetAmphilPartPos<Parameter::REAC_PART,Parameter::TAIL_PART>(up,n_v,param,phil_idx,oil_idx,sDPD,false);
  }
}

void Initializer::MakeHalfSheet(int elem,int& phil_idx,int& oil_idx,double x_leng, double y_leng, double z_leng, const Parameter &param,dpdsystem &sDPD,RNG& rng,bool flag)
{
  const int philic_end_idx = phil_idx + elem;
  const double sign = 2.0 * flag - 1.0;
  const double3 n_v(0.0,sign,0.0);
  const double offset = 0.0001;
  const double area_p_lip = param.L.x*param.L.z / (elem / Parameter::REAC_PART);
  double gr_leng = sqrt(area_p_lip);
  const int grid_n[2] ={static_cast<int>(x_leng / gr_leng),
			static_cast<int>(z_leng / gr_leng)};
  gr_leng = x_leng / grid_n[0];
  int idx[2] = {0};
  bool is_filled = false;
  
  while(phil_idx < philic_end_idx){
    const double3 base(gr_leng*(idx[0] + 0.5) + is_filled*offset,
		       y_leng - sign*Parameter::ALL_UNIT_N*param.b_leng - sign*0.01,
		       gr_leng*(idx[1] + 0.5));
    
    SetAmphilPartPos<0,Parameter::REAC_PART>(base,n_v,param,phil_idx,oil_idx,sDPD,true);
    SetAmphilPartPos<Parameter::REAC_PART,Parameter::TAIL_PART>(base,n_v,param,phil_idx,oil_idx,sDPD,false);
    
    idx[0]++;
    if(idx[0] == grid_n[0]){
      idx[0]=0;
      idx[1]++;
    }
    if(idx[1] == grid_n[1]){
      is_filled = true;
      idx[0]=0;
      idx[1]=0;
    }
  }
}

double Initializer::MakeCylindSheet(int& phil_idx, int& oil_idx, double cylrad, double cyl_Ly, double offset, int cylelem,  const Parameter &param, dpdsystem &sDPD, RNG& rng)
{
  const double lat_leng = sqrt(param.lip_area);

  const int y_elem      = static_cast<int>(cyl_Ly / lat_leng);
  const int disk_elem   = cylelem / y_elem;
  const int res_elem    = cylelem - y_elem * disk_elem;

  const int in_elem     = static_cast<int>(disk_elem * (cylrad - Parameter::ALL_UNIT_N * Parameter::b_leng * 0.5) / (2.0 * cylrad));
  const int out_elem    = disk_elem - in_elem;

  const double d_the_in = 2.0 * M_PI / in_elem;
  const double d_the_out= 2.0 * M_PI / out_elem;

  const double in_rad	= cylrad - Parameter::ALL_UNIT_N * param.b_leng - 0.001;
  const double out_rad	= cylrad + Parameter::ALL_UNIT_N * param.b_leng + 0.001;

  for(int iy = 0; iy < y_elem; iy++){
    const double3 cent(param.hL.x, iy * lat_leng + offset, param.hL.z);
    MakeArcLine(in_elem , d_the_in , cent, in_rad ,  1.0, phil_idx, oil_idx, sDPD, param);
    MakeArcLine(out_elem, d_the_out, cent, out_rad, -1.0, phil_idx, oil_idx, sDPD, param);
  }

  const double3 cent(param.hL.x, y_elem * lat_leng + offset, param.hL.z);
  MakeArcLine(res_elem, d_the_out, cent, out_rad, -1.0, phil_idx, oil_idx, sDPD, param);

  const double last_cylLy = (y_elem - 1) * lat_leng + offset;
  return last_cylLy;
}

void Initializer::MakeArcLine(int elem, double d_the, double3 center, double rad, double sign, int& phil_idx, int& oil_idx, dpdsystem &sDPD, const Parameter &param){
  for(int i=0; i<elem; i++){
    const double cur_the = i * d_the;
    const double3 cent2arc(center.x + rad * cos(cur_the),
			   center.y,
			   center.z + rad * sin(cur_the));
    const double3 n_v(sign * cos(cur_the),
		      0.0,
		      sign * sin(cur_the));
    SetAmphilPartPos<0,Parameter::REAC_PART>(cent2arc, n_v, param, phil_idx, oil_idx, sDPD, true);
    SetAmphilPartPos<Parameter::REAC_PART,Parameter::TAIL_PART>(cent2arc, n_v, param, phil_idx, oil_idx, sDPD, false);
  }
}

//NOTE: rad is the radius of membrane's neutral surface.
bool Initializer::MakeSphLine(const double the, const double3& center, const double rad, const double sign, const double lip_len,
			      const int end_phil_idx, int& phil_idx, int& oil_idx, dpdsystem& sDPD, const Parameter& param)
{
  const double half_thick = 0.5 * Parameter::ALL_UNIT_N * Parameter::b_leng ;
  const double proj_rad = rad * std::sin(the);
  
  double d_phi = lip_len / proj_rad;
  int phi_elem = static_cast<int>(2.0 * M_PI / d_phi);
  d_phi = 2.0 * M_PI / phi_elem;

  for(int i=0; i<phi_elem; i++){
    const double phi = d_phi * i;
    const double3 n_v(sign * std::sin(the) * std::cos(phi), 
		      sign * std::sin(the) * std::sin(phi),
		      sign * std::cos(the));
    const double3 cent2arc = sign * n_v * (rad - sign * half_thick) + center;
    SetAmphilPartPos<0,Parameter::REAC_PART>(cent2arc, n_v, param, phil_idx, oil_idx, sDPD, true);
    SetAmphilPartPos<Parameter::REAC_PART,Parameter::TAIL_PART>(cent2arc, n_v, param, phil_idx, oil_idx, sDPD, false);
    if(phil_idx >= end_phil_idx) return false;
  }

  return true;
}

void Initializer::MakeForEachTheta(const int elem, const double d_the, const double offset, const double rad, const double sign, 
				   const double lip_len, const int end_phil_idx, int& phil_idx, int& oil_idx, dpdsystem& sDPD, 
				   const Parameter& param)
{
  const double3 cent = param.hL;
  for(int i=0; i<elem; i++){
    const double the = d_the * (i + offset);
    const bool flag = MakeSphLine(the, cent, rad, sign, lip_len, end_phil_idx, phil_idx, oil_idx, sDPD, param);
    if(!flag){
      return;
    }
  }
}

void Initializer::MakeSphereSheet(int elem, int& phil_idx, int&oil_idx, const Parameter& param, dpdsystem& sDPD){
  const int end_phil_idx = elem + phil_idx;
  const double lip_len    = std::sqrt(Parameter::lip_area);
  const double half_thick = 0.5 * Parameter::ALL_UNIT_N * Parameter::b_leng;
  const double out_rad = param.binfo.bind_radius + half_thick;
  const double in_rad  = param.binfo.bind_radius - half_thick;
  
  double d_the_out = lip_len / out_rad; 
  double d_the_in  = lip_len / in_rad;  
  
  const int out_the_elem = static_cast<int>(M_PI / d_the_out);
  const int in_the_elem  = static_cast<int>(M_PI / d_the_in);
  
  d_the_out = M_PI / out_the_elem;
  d_the_in  = M_PI / in_the_elem;
  
  //in
  MakeForEachTheta(in_the_elem, d_the_in, 0.0, in_rad, 1.0, lip_len, end_phil_idx, phil_idx, oil_idx, sDPD, param);
  //out
  MakeForEachTheta(out_the_elem, d_the_out, 0.0, out_rad, -1.0, lip_len, end_phil_idx, phil_idx, oil_idx, sDPD, param);

  const int residue = end_phil_idx - phil_idx;
  
  std::cout << "spherical component is " << phil_idx << "." << std::endl;
  std::cout << "residue is " << residue << "." << std::endl;

  if(residue != 0){
    //residue
    MakeForEachTheta(in_the_elem, d_the_in, 0.5, in_rad, 1.0, lip_len, end_phil_idx, phil_idx, oil_idx, sDPD, param);
  }
}

void Initializer::MakeWater(int elem, int& w_idx, const double3& reg_dw, const double3& reg_up, const Parameter &param, dpdsystem &sDPD, RNG& rng, bool flag)
{
  const int w_end_idx = elem + w_idx;
  while(w_idx < w_end_idx){
    const double3 r(rng.Uniform(0) * param.L.x,
		    rng.Uniform(0) * param.L.y,
		    rng.Uniform(0) * param.L.z);
    if(!(((r < reg_up) && (r > reg_dw)) ^ flag)){
      sDPD.SetPosition(w_idx,r);
      sDPD.SetProp(w_idx,Water);
      sDPD.SetLipidUnit(w_idx,-1);
      w_idx++;
    }
  }
}

void Initializer::MakeWaterAsym(int elem, int& w_idx, const double3& reg_dw, const double3& reg_up, const Parameter &param, dpdsystem &sDPD, RNG& rng, const double ratio)
{
  const int up_elem = static_cast<int>(elem * ratio);
  const int dw_elem = elem - up_elem;
  
  const int up_end_idx = elem + up_elem;
  const int dw_end_idx = up_end_idx + dw_elem;
  
  while(w_idx < up_end_idx){
    const double3 r(rng.Uniform(0) * param.L.x,
		    rng.Uniform(0) * param.L.y,
		    rng.Uniform(0) * param.L.z);
    if(!(r < reg_up)){
      sDPD.SetPosition(w_idx,r);
      sDPD.SetProp(w_idx,Water);
      sDPD.SetLipidUnit(w_idx,-1);
      w_idx++;
    }
  }
  
  while(w_idx < dw_end_idx){
    const double3 r(rng.Uniform(0) * param.L.x,
		    rng.Uniform(0) * param.L.y,
		    rng.Uniform(0) * param.L.z);
    if(!(r > reg_dw)){
      sDPD.SetPosition(w_idx,r);
      sDPD.SetProp(w_idx,Water);
      sDPD.SetLipidUnit(w_idx,-1);
      w_idx++;
    }
  }
}

void Initializer::MakeWaterCyl(int elem, int& w_idx, double cyl_rad, double3 cent, const double* y_reg, const Parameter &param, dpdsystem &sDPD, RNG& rng)
{
  const double out_rad = cyl_rad + Parameter::ALL_UNIT_N * param.b_leng;
  const double in_rad  = cyl_rad - Parameter::ALL_UNIT_N * param.b_leng;
  const int w_end_idx = elem + w_idx;
  while(w_idx < w_end_idx){
    const double3 r(rng.Uniform(0) * param.L.x,
		    rng.Uniform(0) * param.L.y,
		    rng.Uniform(0) * param.L.z);
    const double3 cent2r = r - cent;
    const double dist = sqrt(cent2r.x * cent2r.x + cent2r.z * cent2r.z);
    if((dist > out_rad) || (dist < in_rad) || (r.y > y_reg[1]) || (r.y < y_reg[0])){
      sDPD.SetPosition(w_idx, r);
      sDPD.SetProp(w_idx, Water);
      sDPD.SetLipidUnit(w_idx, -1);
      w_idx++;
    }
  }
}

void Initializer::MakeWaterSphere(const int elem, int& w_idx, const int in_elem, const Parameter& param, dpdsystem& sDPD, RNG& rng){
  const int w_end_idx = elem + w_idx;
  const double half_thick = Parameter::ALL_UNIT_N * Parameter::b_leng;
  const int out_elem = elem - in_elem;

  int in_count = 0, out_count = 0;
  
  while(w_idx < w_end_idx){
    const double3 r(rng.Uniform(0) * param.L.x,
		    rng.Uniform(0) * param.L.y,
		    rng.Uniform(0) * param.L.z);
    const double3 cent2r = r - param.hL;
    const double dist = cent2r.norm2();
    if( dist - param.binfo.bind_radius > half_thick){
      if(out_count < out_elem){
	sDPD.SetPosition(w_idx, r);
	sDPD.SetProp(w_idx, Water);
	sDPD.SetLipidUnit(w_idx, -1);
	w_idx++;
	out_count++;
      }
    }else if (dist - param.binfo.bind_radius < -half_thick){
      if(in_count < in_elem){
	sDPD.SetPosition(w_idx, r);
	sDPD.SetProp(w_idx, Water);
	sDPD.SetLipidUnit(w_idx, -1);
	w_idx++;
	in_count++;
      }
    }
  }
}

void Initializer::MakeEmbedTail(int elem, int& phil_idx, int& oil_idx, const double3& reg_dw, const double3& reg_up, const Parameter &param,dpdsystem &sDPD,RNG& rng)
{
  const int oil_end_idx = oil_idx + elem;
  const double3 box_leng = reg_up - reg_dw;
  
#if 1
  const double offset = 0.0001, area_p_lip = param.L.x * param.L.z / (elem / Parameter::TAIL_PART);
  double gr_leng = sqrt(area_p_lip);
  const int grid_n[2] ={static_cast<int>(box_leng.x / gr_leng),
			static_cast<int>(box_leng.z / gr_leng)};
  gr_leng = box_leng.x / grid_n[0];
  int idx[2] = {0}; bool is_filled = false;
  
  while(oil_idx < oil_end_idx) {
    double3 base(gr_leng * (idx[0] + 0.5) + is_filled * offset,
		 0.5 * Parameter::TAIL_PART * Parameter::b_leng,
		 gr_leng * (idx[1] + 0.5));
    base += reg_dw;
    const double3 n_v(0.0, 1.0, 0.0);
    SetAmphilPartPos<Parameter::REAC_PART,Parameter::TAIL_PART>(base, n_v, param, phil_idx, oil_idx, sDPD, false);
    idx[0]++;
    if(idx[0] == grid_n[0]){
      idx[0]=0;
      idx[1]++;
    }
    if(idx[1] == grid_n[1]){
      is_filled = true;
      idx[0]=0;
      idx[1]=0;
    }
  }
  
#else

  while(oil_idx < oil_end_idx){
    double3 r(box_leng.x * rng.Uniform(0),
	      0.0,
	      box_leng.z * rng.Uniform(0));
    r += reg_dw;
      
    const double3 n_v(0.0,1.0,0.0);
    SetAmphilPartPos<Parameter::REAC_PART,Parameter::TAIL_PART>(r,n_v,param,phil_idx,oil_idx,sDPD,false);
  }

#endif
}

void Initializer::MakeHeadInWater(int elem, int& phil_idx, int& oil_idx, const double3& reg_dw, const double3& reg_up, const Parameter &param, dpdsystem &sDPD,RNG& rng)
{
  const int phil_end_idx = elem + phil_idx;
  const double3 box_leng = reg_up - reg_dw;
    
  while(phil_idx < phil_end_idx){
    double3 r(rng.Uniform(0) * box_leng.x,
	      rng.Uniform(0) * box_leng.y,
	      rng.Uniform(0) * box_leng.z);
    r += reg_dw;
      
    double3 n_v(rng.Uniform(0),
		rng.Uniform(0),
		rng.Uniform(0));
    n_v /= n_v.norm2();
    SetAmphilPartPos<0,Parameter::REAC_PART>(r,n_v,param,phil_idx,oil_idx,sDPD,true);
  }
}

void Initializer::MakeHeadInWaterCyl(int elem, int& phil_idx, int& oil_idx, double cyl_rad, double3 cent, double z_hei, bool flag, const Parameter &param, dpdsystem &sDPD, RNG& rng)
{
  const int phil_end_idx = elem + phil_idx;

  while(phil_idx < phil_end_idx){
    double3 r(rng.Uniform(0) * param.L.x,
	      rng.Uniform(0) * param.L.y,
	      rng.Uniform(0) * param.L.z);
    const double3 cent2r = r - cent;
    const double dist = sqrt(cent2r.x*cent2r.x + cent2r.y*cent2r.y);
    if(r.z < z_hei){
      if( (dist < cyl_rad) ^ flag){
	double3 n_v(rng.Uniform(0),
		    rng.Uniform(0),
		    rng.Uniform(0));
	n_v /= n_v.norm2();
	SetAmphilPartPos<0,Parameter::REAC_PART>(r, n_v, param, phil_idx, oil_idx, sDPD, true);
      }
    }
  }
}

//NOTE: Head molecules are disperted in vesicle.
void Initializer::MakeHeadInWaterSph(int elem, int& phil_idx, int& oil_idx, const Parameter &param, dpdsystem &sDPD, RNG& rng)
{
  const int phil_end_idx = elem + phil_idx;
  const double half_thick = Parameter::ALL_UNIT_N * Parameter::b_leng;
  
  while(phil_idx < phil_end_idx){
    const double3 r(rng.Uniform(0) * param.L.x,
		    rng.Uniform(0) * param.L.y,
		    rng.Uniform(0) * param.L.z);
    const double3 cent2r = r - param.hL;
    const double dist = cent2r.norm2();
    if(dist - param.binfo.bind_radius < -half_thick){
      double3 n_v(rng.Uniform(0),
		  rng.Uniform(0),
		  rng.Uniform(0));
      n_v /= n_v.norm2();
      SetAmphilPartPos<0,Parameter::REAC_PART>(r, n_v, param, phil_idx, oil_idx, sDPD, true);
    }
  }
}

void Initializer::MakeOilbulk(int elem, int& phil_idx, int& oil_idx, const double3& reg_dw, const double3& reg_up, const Parameter& param, dpdsystem& sDPD, RNG& rng){
  const int o_end_idx = elem + oil_idx;
  const double3 box_leng = reg_up - reg_dw;
  while(oil_idx < o_end_idx){
    double3 r(rng.Uniform(0) * box_leng.x,
	      rng.Uniform(0) * box_leng.y,
	      rng.Uniform(0) * box_leng.z);
    r += reg_dw;

    double3 n_v(rng.Uniform(0),
		rng.Uniform(0),
		rng.Uniform(0));
    n_v /= n_v.norm2();
    SetAmphilPartPos<Parameter::REAC_PART,Parameter::TAIL_PART>(r,n_v,param,phil_idx,oil_idx,sDPD,false);
  }
}

void Initializer::GenDiskLam(dpdsystem &sDPD,const Parameter& param,RNG& rng,int mode){
  int phil_idx=0, oil_idx=param.hN, w_idx=param.hN+param.bN;
    
  //amphiphiles
  const double3 offset(1.5, 0.0, 1.5);
  const double z_strp_leng = 5.0;
  MakeDisk(param.ampN * Parameter::HYPHIL_N, phil_idx, oil_idx, offset, param, sDPD);
  if(mode == NORMAL){
    assert(param.tailN == param.headN);
    assert(phil_idx == param.hN); assert(oil_idx == param.hN + param.bN);      
  }else if(mode == TAIL_EMBEDDED){
    assert(param.headN < param.tailN);
    const int elem = (Parameter::TAIL_PART > Parameter::HYPHOB_N) ? (param.tailN - param.headN)*Parameter::HYPHOB_N : (param.tailN - param.headN)*Parameter::TAIL_PART;
    const double3 reg_dw(offset.x, param.hL.y - Parameter::ALL_UNIT_N * param.b_leng, offset.z);
    const double3 reg_up(param.L.x - offset.x, param.hL.z + Parameter::ALL_UNIT_N * param.b_leng, offset.z + z_strp_leng);
    MakeEmbedTail(elem, phil_idx, oil_idx, reg_dw, reg_up, param, sDPD,rng);
  }else if(mode == HEAD_DISPERT){
    assert(param.headN > param.tailN);
    const int elem = (Parameter::REAC_PART > Parameter::HYPHIL_N) ? (param.headN - param.tailN)*Parameter::HYPHIL_N : (param.headN - param.tailN)*Parameter::REAC_PART;
    const double3 reg_dw(0.0);
    const double3 reg_up(param.L.x, param.L.y, param.L.z);
    MakeHeadInWater(elem, phil_idx, oil_idx, reg_dw, reg_up, param, sDPD,rng);
  }

  //water
  const double3 reg_dw(offset.x ,param.hL.y - Parameter::ALL_UNIT_N*param.b_leng - 0.01,offset.z );
  const double3 reg_up(param.L.x,param.hL.y + Parameter::ALL_UNIT_N*param.b_leng + 0.01,param.L.z);
  MakeWater(param.wN,w_idx,reg_dw,reg_up,param,sDPD,rng,false);
  assert(w_idx == Parameter::SYS_SIZE);
}

void Initializer::GenPlaneLam(dpdsystem &sDPD,const Parameter& param,RNG& rng,int mode){
  int phil_idx=0, oil_idx=param.hN, w_idx=param.hN+param.bN;

  //amphiphiles
  MakeSheet(param.ampN*Parameter::HYPHIL_N, phil_idx, oil_idx, param.L.x, param.hL.y, param.L.z, param, sDPD, rng);
  if(mode == TAIL_EMBEDDED){
    assert(param.headN < param.tailN);
    const int elem = (Parameter::TAIL_PART > Parameter::HYPHOB_N) ? (param.tailN - param.headN)*Parameter::HYPHOB_N : (param.tailN - param.headN)*Parameter::TAIL_PART;
    const double3 reg_dw(0.0,param.hL.y - Parameter::ALL_UNIT_N*param.b_leng,0.0);
    const double3 reg_up(param.L.x,param.hL.y + Parameter::ALL_UNIT_N*param.b_leng,param.L.z);
    MakeEmbedTail(elem, phil_idx, oil_idx, reg_dw, reg_up, param, sDPD,rng);
    assert(phil_idx == param.hN); assert(oil_idx == param.hN + param.bN);      
  }else if(mode == NORMAL){
    assert(param.headN == param.tailN);
    assert(phil_idx == param.hN); assert(oil_idx == param.hN + param.bN);      
  }else if(mode == HEAD_DISPERT){
    assert(param.headN > param.tailN);
    const float ratio = 0.8; //upper region : downer region = 8 : 2
    const int up_h_elem = static_cast<int>((param.headN - param.tailN) * ratio);
    const int dw_h_elem = (param.headN - param.tailN) - up_h_elem; 
    const int up_elem = (Parameter::REAC_PART > Parameter::HYPHIL_N) ? up_h_elem*Parameter::HYPHIL_N : up_h_elem*Parameter::REAC_PART;
    const int dw_elem = (Parameter::REAC_PART > Parameter::HYPHIL_N) ? dw_h_elem*Parameter::HYPHIL_N : dw_h_elem*Parameter::REAC_PART;
    
    double3 reg_up[2], reg_dw[2];
    reg_up[0].x = reg_up[0].z = 0.0; 
    reg_up[0].y = param.hL.y + Parameter::ALL_UNIT_N*param.b_leng + 0.01;
    reg_up[1]   = param.L; reg_up[1].y -= Parameter::ALL_UNIT_N*param.b_leng;
    reg_dw[0].x = reg_dw[0].y = reg_dw[0].z = 0.0; reg_dw[0].y += Parameter::ALL_UNIT_N*param.b_leng;
    reg_dw[1].x = param.L.x; 
    reg_dw[1].y = param.hL.y - Parameter::ALL_UNIT_N*param.b_leng - 0.01; 
    reg_dw[1].z = param.L.z;
    
    MakeHeadInWater(up_elem, phil_idx, oil_idx, reg_up[0], reg_up[1], param, sDPD, rng);
    MakeHeadInWater(dw_elem, phil_idx, oil_idx, reg_dw[0], reg_dw[1], param, sDPD, rng);
  }

  //water
  const double3 reg_dw(0.0, param.hL.y - Parameter::ALL_UNIT_N*param.b_leng - 0.01, 0.0);
  const double3 reg_up(param.L.x, param.hL.y + Parameter::ALL_UNIT_N*param.b_leng + 0.01, param.L.z);

  //sym
  MakeWater(param.wN, w_idx, reg_dw, reg_up, param, sDPD, rng, false);
  
  //asym
  //const double ratio = 0.6;
  //MakeWaterAsym(param.wN, w_idx, reg_dw, reg_up, param, sDPD, rng, ratio);
  assert(w_idx == Parameter::SYS_SIZE);
}

void Initializer::GenAsynPlaneLam(dpdsystem &sDPD,const Parameter& param,RNG& rng,int mode){
  int phil_idx = 0, oil_idx = param.hN, w_idx = param.hN+param.bN;
  const int tens_zero_N = static_cast<int>(param.L.x * param.L.z / param.lip_area);
  const int res_N = param.ampN - tens_zero_N;
  assert(res_N > 0);

  MakeHalfSheet(tens_zero_N * Parameter::HYPHIL_N, phil_idx, oil_idx, param.L.x, param.hL.y, param.L.z, param, sDPD, rng, true);
  MakeHalfSheet(res_N * Parameter::HYPHIL_N, phil_idx, oil_idx, param.L.x, param.hL.y, param.L.z, param, sDPD, rng, false);

  if(mode == NORMAL) {
    assert(phil_idx == param.hN); assert(oil_idx == param.hN + param.bN);
  } else if(mode == TAIL_EMBEDDED) {
    assert(param.headN < param.tailN);
    const int elem = (Parameter::TAIL_PART > Parameter::HYPHOB_N) ? (param.tailN - param.headN) * Parameter::HYPHOB_N : (param.tailN - param.headN) * Parameter::TAIL_PART;
    const double3 reg_dw(0.0, param.hL.y - Parameter::ALL_UNIT_N * param.b_leng, 0.0);
    const double3 reg_up(param.L.x,param.hL.y + Parameter::ALL_UNIT_N*param.b_leng,param.L.z);
    MakeEmbedTail(elem, phil_idx, oil_idx, reg_dw, reg_up, param, sDPD, rng);
    assert(phil_idx == param.hN); assert(oil_idx == param.hN + param.bN);
  } else {
    std::cerr << "Unknown mode \n";
    std::exit(1);
  }
  
  //water
  const double3 reg_dw(0.0      ,param.hL.y - Parameter::ALL_UNIT_N*param.b_leng - 0.01      ,0.0);
  const double3 reg_up(param.L.x,param.hL.z + Parameter::ALL_UNIT_N*param.b_leng + 0.01,param.L.z);
  MakeWater(param.wN,w_idx,reg_dw,reg_up,param,sDPD,rng,false);
  assert(w_idx == Parameter::SYS_SIZE);
}

void Initializer::GenStrip(dpdsystem &sDPD,const Parameter& param, RNG& rng,int mode){
  assert(mode < 2);
  int phil_idx=0, oil_idx=param.hN, w_idx=param.hN+param.bN;

  //amphiphiles
  const double strL = 0.7 * param.L.x;
  MakeSheet(param.ampN*Parameter::HYPHIL_N, phil_idx, oil_idx, strL, param.hL.y, param.L.z, param, sDPD,rng);
  if(mode == TAIL_EMBEDDED){
    assert(param.headN < param.tailN);
    const int elem = (Parameter::TAIL_PART > Parameter::HYPHOB_N) ? (param.tailN - param.headN)*Parameter::HYPHOB_N : (param.tailN - param.headN)*Parameter::TAIL_PART;
    const double3 reg_dw(0.0,param.hL.y - Parameter::ALL_UNIT_N*param.b_leng,0.0);
    const double3 reg_up(strL,param.hL.y + Parameter::ALL_UNIT_N*param.b_leng,param.L.z);
    MakeEmbedTail(elem, phil_idx, oil_idx, reg_dw, reg_up, param, sDPD,rng);
    assert(phil_idx == param.hN); assert(oil_idx == param.hN + param.bN);      
  }else if(mode == NORMAL){
    assert(param.tailN == param.headN);
    assert(phil_idx == param.hN); assert(oil_idx == param.hN + param.bN);
  }
    
  //water
  const double3 reg_dw(0.0 ,param.hL.y - Parameter::ALL_UNIT_N*param.b_leng - 0.01      ,0.0);
  const double3 reg_up(strL,param.hL.y + Parameter::ALL_UNIT_N*param.b_leng + 0.01,param.L.z);
  MakeWater(param.wN,w_idx,reg_dw,reg_up,param,sDPD,rng,false);
  assert(w_idx == Parameter::SYS_SIZE);
}

void Initializer::GenRandomSolute(dpdsystem &sDPD,const Parameter& param,RNG& rng){
  assert(param.tailN == param.headN);
  int phil_idx=0, oil_idx=param.hN, w_idx=param.hN+param.bN;
    
  //amphiphiles
  MakeRandomDispert(param.ampN*Parameter::HYPHIL_N,phil_idx,oil_idx,param,sDPD,rng);
  assert(phil_idx == param.hN); assert(oil_idx == param.hN + param.bN);
    
  //waters
  const double3 zero(0.0);
  MakeWater(param.wN,w_idx,zero,zero,param,sDPD,rng,false);
  assert(w_idx == Parameter::SYS_SIZE);
}  

void Initializer::GenOilDrop(dpdsystem &sDPD,const Parameter& param,RNG& rng){
  assert(param.tailN <= param.headN);

  int phil_idx=0, oil_idx=param.hN, w_idx=param.hN+param.bN;
  const double	dropR	     = Pow_1_n(param.bN*0.25/M_PI,3);

  int tail_idx=0;
  while(tail_idx < param.tailN){
    double3 n_v(rng.Uniform(0),rng.Uniform(0),rng.Uniform(0));
    n_v /= n_v.norm2();
    const double3 r(rng.Uniform(0) * param.L.x,
		    rng.Uniform(0) * param.L.y,
		    rng.Uniform(0) * param.L.z);
    const double3 cent2r = r - param.hL;
    
    if(cent2r.norm2() < dropR - Parameter::ALL_UNIT_N * param.b_leng){
      SetAmphilPartPos<Parameter::REAC_PART,Parameter::TAIL_PART>(r,n_v,param,phil_idx,oil_idx,sDPD,false);
      tail_idx++;
    }
  }

  int head_idx=0;
  while(head_idx < param.headN){
    double3 n_v(rng.Uniform(0),rng.Uniform(0),rng.Uniform(0));
    n_v /= n_v.norm2();
    const double3 r(rng.Uniform(0) * param.L.x,
		    rng.Uniform(0) * param.L.y,
		    rng.Uniform(0) * param.L.z);
    const double3 cent2r = r - param.hL;

    if(cent2r.norm2() >= dropR - Parameter::ALL_UNIT_N * param.b_leng){
      SetAmphilPartPos<0,Parameter::REAC_PART>(r,n_v,param,phil_idx,oil_idx,sDPD,true);
      head_idx++;
    }
  }

  const int	w_end_idx    = Parameter::SYS_SIZE;
  while(w_idx < w_end_idx){
    const double3 r(rng.Uniform(0) * param.L.x,
		    rng.Uniform(0) * param.L.y,
		    rng.Uniform(0) * param.L.z);
    const double3 cent2r = r - param.hL;
    
    if(cent2r.norm2() >= dropR - Parameter::ALL_UNIT_N * param.b_leng){
      sDPD.SetPosition(w_idx,r);
      sDPD.SetProp(w_idx,Water);
      sDPD.SetLipidUnit(w_idx,-1);
      w_idx++;
    }
  }
  
  assert(phil_idx == param.hN); assert(oil_idx == param.hN+param.bN);
}

void Initializer::GenOilCuboid(dpdsystem &sDPD, const Parameter& param, RNG& rng){
  int phil_idx=0, oil_idx=param.hN, w_idx=param.hN+param.bN;
  const double offset = Parameter::ALL_UNIT_N * param.b_leng;
  const double oil_height = static_cast<double>(param.bN) / static_cast<double>(Parameter::SYS_SIZE) * param.L.y + offset;
  
  const double3 origin(0.0);
  const double3 reg_dw(0.0,offset,0.0);
  const double3 reg_up(param.L.x,oil_height,param.L.z);

  const double3 h_reg_dw(0.0,oil_height,0.0);
  const double3 h_reg_up(param.L.x, param.L.y-offset, param.L.z);
  
  const int o_elem = (Parameter::TAIL_PART > Parameter::HYPHOB_N) ? param.tailN * Parameter::HYPHOB_N : param.tailN * Parameter::TAIL_PART;
  const int h_elem = (Parameter::REAC_PART > Parameter::HYPHIL_N) ? param.headN * Parameter::HYPHIL_N : param.headN * Parameter::REAC_PART;

  MakeOilbulk(o_elem,phil_idx,oil_idx,reg_dw,reg_up,param,sDPD,rng);
  MakeHeadInWater(h_elem,phil_idx,oil_idx,h_reg_dw,h_reg_up,param,sDPD,rng);
  MakeWater(param.wN,w_idx,origin,reg_up,param,sDPD,rng,false);

  assert(phil_idx == param.hN);
  assert(oil_idx == param.hN + param.bN);
  assert(w_idx == Parameter::SYS_SIZE);
}

void Initializer::GenCylindLam(dpdsystem &sDPD, const Parameter& param, RNG& rng, int mode){
  int phil_idx = 0, oil_idx = param.hN, w_idx = param.hN + param.bN;
  const double cyl_Ly = param.L.y * 0.85;
  const int cyl_elem  = param.ampN;
  
  //NOTE: This value shold be tuned depending on cyl_Lz.
  const double cyl_R  = cyl_elem * param.lip_area * 0.25 / M_PI / cyl_Ly + 3.0;

  std::cout << "The radius of bilayer cylinder is " << cyl_R << std::endl;

  const double3 cent(param.hL.x, 0.0, param.hL.z);

  assert(2.0 * (cyl_R + Parameter::ALL_UNIT_N * Parameter::b_leng) < param.L.x);
  assert(cyl_R > Parameter::ALL_UNIT_N * Parameter::b_leng);

  //amphiphiles
  const double last_cyl_Ly = MakeCylindSheet(phil_idx, oil_idx, cyl_R, cyl_elem, cyl_Ly, 0.0, param, sDPD, rng);
  if(mode == NORMAL){
    assert(param.headN == param.tailN);
    assert(phil_idx == param.hN      ); assert(oil_idx == param.hN + param.bN);
  }else if(mode == HEAD_DISPERT){
    assert(param.headN > param.tailN);
    const int elem		= (Parameter::REAC_PART > Parameter::HYPHIL_N) ? (param.headN - param.tailN) * Parameter::HYPHIL_N : (param.headN - param.tailN) * Parameter::REAC_PART;
    const double ratio		= 0.8;
    const double in_rad		= cyl_R - Parameter::ALL_UNIT_N * param.b_leng;
    const int    in_elem	= static_cast<int>(elem * ratio);
    const double out_rad	= cyl_R + Parameter::ALL_UNIT_N * param.b_leng;
    const int    out_elem	= elem - in_elem;
    
    MakeHeadInWaterCyl(in_elem , phil_idx, oil_idx, in_rad , cent, cyl_Ly, false, param, sDPD, rng);
    MakeHeadInWaterCyl(out_elem, phil_idx, oil_idx, out_rad, cent, cyl_Ly, true , param, sDPD, rng);
    assert(phil_idx == param.hN      ); assert(oil_idx == param.hN + param.bN);    
  }
  
  //water
  const double y_reg[2] = {0.0, last_cyl_Ly};
  MakeWaterCyl(param.wN, w_idx, cyl_R, cent, y_reg, param, sDPD, rng);
  assert(w_idx == Parameter::SYS_SIZE);
}

void Initializer::GenSphere(dpdsystem& sDPD, const Parameter& param, RNG& rng){
  int phil_idx = 0, oil_idx = param.hN, w_idx = param.hN + param.bN;
  
  assert(param.headN >= param.tailN);
  
  // const double red_vol = 0.4;
  // const double all_area = Parameter::lip_area * param.ampN * 0.5;
  // const double vol = red_vol * 4.0 / 3.0 * M_PI * (all_area * 0.25 / M_PI) * sqrt(all_area * 0.25 / M_PI);
  // const int in_elem = static_cast<int>(vol * 3.0);
  const int in_elem = 10000;

  assert(in_elem < param.wN);

  MakeSphereSheet(Parameter::HYPHIL_N * param.ampN, phil_idx, oil_idx, param, sDPD);
  MakeWaterSphere(param.wN, w_idx, in_elem, param, sDPD, rng);

  if(param.headN > param.tailN){
    const int disp_elem = (Parameter::REAC_PART > Parameter::HYPHIL_N) ? (param.headN - param.tailN) * Parameter::HYPHIL_N : (param.headN - param.tailN) * Parameter::REAC_PART;
    MakeHeadInWaterSph(disp_elem, phil_idx, oil_idx, param, sDPD, rng);
  }
  
  assert(phil_idx == param.hN);      
  assert(oil_idx == param.hN + param.bN);
  assert(w_idx == Parameter::SYS_SIZE);
}

void Initializer::LoadRestartConfig(dpdsystem& sDPD, ChemInfo& cheminfo, const Parameter& param, std::ifstream& fin) const {
  int id = 0;
  int w_n = 0, h_n = 0, b_n = 0;
  double3 buf_p(0.0), buf_v(0.0);
  bool buf_p_chem = false, buf_l_chem = false;
  int buf_l_unit = -1, buf_l_idx = -1, buf_p_idx = -1, buf_prop = -1;
  
  while(true){
    fin >> buf_p.x >> buf_p.y >> buf_p.z >> buf_v.x >> buf_v.y >> buf_v.z >> buf_prop >> buf_p_chem
	>> buf_l_chem >> buf_l_unit >> buf_l_idx >> buf_p_idx;
    if(fin.eof() ) break;
    
    sDPD.pr[id]			= buf_p;
    sDPD.pv[id]			= buf_v;
    sDPD.prop[id]		= static_cast<par_prop>(buf_prop);
    cheminfo.prtcl_chem[id]	= buf_p_chem;
    cheminfo.lipid_idx[id]	= buf_l_idx;
    cheminfo.lipid_unit[id]	= buf_l_unit;
    cheminfo.part_idx[id]	= buf_p_idx;
    
    id++;

    if(buf_prop == Water){
      w_n++;
    }else if(buf_prop == Hyphil){
      h_n++;
    }else if(buf_prop == Hyphob){
      b_n++;
    }
    
    if(id > Parameter::SYS_SIZE){
      std::cerr << "The number of particle is larger than SYS_SIZE \n";
      std::exit(1);
    }
  }
    
  if(id != Parameter::SYS_SIZE){
    std::cerr << "The number of particle is smaller than SYS_SIZE \n";
    std::exit(1);
  }

  NumberIsMatched(w_n, param.wN, "Water");
  NumberIsMatched(h_n, param.hN, "Hyphil");
  NumberIsMatched(b_n, param.bN, "Hyhpob");
}

void Initializer::NumberIsMatched(const int r_v, const int l_v, const std::string& mess) const {
  if(r_v != l_v){
    std::cerr << "Unmatched : particle name is " << mess			<< std::endl;;
    std::cerr << "Loaded particle number from init_config.txt is " << r_v	<< std::endl;
    std::cerr << "Loaded particle number from macro_IF.txt is " << l_v		<< std::endl;
    std::exit(1);
  }
}

void Initializer::SetBondedParameter(dpdsystem& sDPD, ChemInfo& cheminfo, const Parameter& param){
  //cheminfo.lip_elem_idx cheminfo.head_elem_idx cheminfo.tail_elem_idx -> ChemManager
  
  for(int i=0; i<param.ampN; i++) cheminfo.lipid_chem[i] = false;
  
  for(int i=0; i<Parameter::SYS_SIZE; i++){
    const int l_unit = cheminfo.lipid_unit[i];
    const int l_idx  = cheminfo.lipid_idx[i];
    if( (l_unit == Parameter::REAC_PART) || (l_unit == Parameter::REAC_PART - 1) ){
      const bool prtcl_chem = cheminfo.prtcl_chem[i];
      
      if(prtcl_chem){
	if(l_idx >= 0 && l_idx < param.ampN){
	  cheminfo.lipid_chem[l_idx] = true;
	}else{
	  std::cerr << "lipid idx is out-of-range. \n";
	  std::cerr << "lipid idx: " << l_idx << std::endl;
	  std::exit(1);
	}
      }
    }
  }
}

void Initializer::CheckRestartConfigIsValid(const dpdsystem& sDPD, const ChemInfo& cheminfo) const {
  for(int i=0; i<Parameter::SYS_SIZE; i++){
    const int l_unit = cheminfo.lipid_unit[i];
    const int l_idx  = cheminfo.lipid_idx[i];
    const bool prtcl_chem = cheminfo.prtcl_chem[i];
    if( (l_unit == Parameter::REAC_PART) || (l_unit == Parameter::REAC_PART - 1) ){
      assert(cheminfo.lipid_chem[l_idx] == prtcl_chem);
    }
  }

  sDPD.CheckChemInfoInited();
  sDPD.CheckParticleInited();
}

void Initializer::GenParticles(dpdsystem &sDPD, ChemManager& chemmanage, ChemInfo& cheminfo, const Parameter& param, RNG& rng)
{
  std::string fname = param.cur_dir + "/init_config.txt";
  std::ifstream fin(fname.c_str() );

  ClearLipidIdx(sDPD);
  
  if (fin) {
    LoadRestartConfig(sDPD, cheminfo, param, fin);
    SetBondedParameter(sDPD, cheminfo, param);
    chemmanage.RegistLipidIdx(param, cheminfo);
    CheckRestartConfigIsValid(sDPD, cheminfo);
  } else {
    std::cerr << "init_config.txt not found. \n";
    std::cerr << "initial configuration is generated. \n";
    const int bond_n = param.ampN;
    //const int bond_n = 0;
    SetLipidId(sDPD, param, bond_n);
    InitVeloc(sDPD, param, rng);
    InitConfig(sDPD, param, rng);
    CheckInited(sDPD, param);
  }
}

void Initializer::CheckInited(const dpdsystem &sDPD, const Parameter& param) const {
  sDPD.CheckChemInfoInited();
  sDPD.CheckParticleInited();
    
  assert(h_array_id   == param.headN * Parameter::REAC_PART);
  assert(t_array_id   == param.tailN * Parameter::TAIL_PART);
  assert(amp_array_id == param.ampN * Parameter::ALL_UNIT_N);
}
