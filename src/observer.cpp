#include "observer.hpp"
#include "bucket_sorter.hpp"
#include "force_calculator.hpp"
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include <algorithm>
#include <omp.h>
#include <numeric>

Observer::Observer(const Parameter& param) {
  loc_tempera.resize(param.grid_numb[1], 0.0);
  loc_dense.resize(param.grid_numb[1], 0.0);
  loc_vel.resize(param.grid_numb[2], double3(0.0));
  tail_cm_pos.resize((param.tailN > param.headN) ? param.tailN : param.headN, double3(0.0));

  // for height calculation
  cut_grid = static_cast<int>(param.L.x / cut_r);
  cut_r = param.L.x / cut_grid;
  i_cut_r = 1.0 / cut_r;
  hei.resize(cut_grid * cut_grid, 0.0);
  hei_elem.resize(cut_grid * cut_grid, 0);
}

Observer::~Observer() {
  for (int i = 0; i < NUM_FILE; i++) fclose(fp[i]);
}

std::string Observer::Type2Fname(const int type, const Parameter& param) {
  switch (type) {
  case MACRO_VAL:
    return param.cur_dir + "/mval.txt";
  case PTCL_CONFIG:
    return param.cur_dir + "/traject.xyz";
  case LOCAL_VAL:
    return param.cur_dir + "/mlocval.txt";
  case PRESSURE:
    return param.cur_dir + "/pressure.txt";
  case GYRATION:
    return param.cur_dir + "/gyration.txt";
  case CONFIG_TEMP:
    return param.cur_dir + "/configtemp.txt";
  case FINAL_CONFIG:
    return param.cur_dir + "/fin_config.xyz";
  case HEIGHT_DIST:
    return param.cur_dir + "/height_dist.txt";
  default:
    std::cerr << "Unknown file type.\n";
    std::exit(1);
  }
}

void Observer::Initialize(const Parameter& param) {
  for (int i = 0; i < NUM_FILE; i++)
    fp[i] = xfopen(Type2Fname(i, param).c_str(), "w");
  fprintf(fp[MACRO_VAL], "#kT\tD\tThick\tS\n");
  fprintf(fp[LOCAL_VAL], "#height\tkBT\tdensity\n");
  fprintf(fp[HEIGHT_DIST], "%d %.10g %.10g %.10g\n", cut_grid, cut_r, param.L.x, param.L.z);
}

double Observer::CalcKinTempera(const dpdsystem &sDPD) {
  double kin_temp = std::accumulate(sDPD.pv, sDPD.pv + Parameter::sys_size, 0.0,
				    [](const double sum, const double3& val) {return sum + val * val;});
  kin_temp /= 3. * Parameter::sys_size;
  return kin_temp;
}

double Observer::CalcDiffs(const dpdsystem& sDPD) {
  double difsum = std::accumulate(sDPD.delta_sumr, sDPD.delta_sumr + Parameter::sys_size, 0.0,
				  [](const double sum, const double3& val) {return sum + val * val;});
  difsum /= Parameter::sys_size;
  return difsum;
}

double Observer::CalcGyration(const dpdsystem& sDPD, const Parameter& param) {
  //search base position
  double3 base_pos(0.0);
  const double3* t_r = sDPD.pr;
  const par_prop* t_prp = sDPD.prop;
  const int base_sample = 4;
  int baseN = 1, bef_hash = -1;
  for(int i=0; i<Parameter::sys_size; i++){
    const int prtcl_idx = sDPD.GetPrtclIdx(i);
    const double3 pos = t_r[prtcl_idx];
    const par_prop prp = t_prp[prtcl_idx];
    
    if (prp != Water) {
      const int temp_hash = B_sorter::GenHash(pos, param);
      if (temp_hash == bef_hash) {
	base_pos += pos;
	baseN++;
	if (baseN == base_sample) break;
      } else {
	baseN = 1;
	base_pos *= 0.0;
	base_pos += pos;
	bef_hash  = temp_hash;
      }
    }
  }
  base_pos *= 1.0 / base_sample;
  
  //calc center of mass
  double3 cm_pos(0.0);
  for(int i = 0; i < Parameter::sys_size; i++) {
    double3 pos = t_r[i];
    const par_prop prp = t_prp[i];
    if (prp != Water) {
      pos -= base_pos;
      F_calculator::MinImage(pos, param);
      cm_pos += pos;
    }
  }
  cm_pos /= (param.hN + param.bN);
  
  //calc gyration tensor
  double gyr[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  for(int i = 0; i < Parameter::sys_size; i++) {
    double3 pos = t_r[i];
    const par_prop prp = t_prp[i];
    if (prp != Water) {
      pos -= base_pos;
      F_calculator::MinImage(pos,param);
      pos -= cm_pos;
      F_calculator::MinImage(pos,param);
      for (int j = 0; j < 3; j++)
	for (int k = j; k < 3; k++)
	  gyr[j][k] += pos[j] * pos[k];
    }
  }
  const double ihbN = 1.0 / (param.hN + param.bN);
  for (int j = 0; j < 3; j++)
    for (int k = j; k < 3; k++)
      gyr[j][k] *= ihbN;
  gyr[1][0] = gyr[0][1];
  gyr[2][0] = gyr[0][2];
  gyr[2][1] = gyr[1][2];

  //diagonalize
  using namespace Eigen;
  Matrix3d gyrmat;
  for (int i =0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      gyrmat(i,j) = gyr[i][j];

  SelfAdjointEigenSolver<Matrix3d> es(gyrmat);
  if (es.info() != Success) std::abort();
  Vector3d tempvec = es.eigenvalues();
  double egvals[3] = {tempvec(0), tempvec(1), tempvec(2)};
  std::sort(egvals, egvals + 3);
  
  const double gyr_rad  = egvals[0] + egvals[1] + egvals[2];
  const double aspheric = 1.5 * egvals[2] - 0.5 * gyr_rad;
  const double acylind  = egvals[1] - egvals[0];
  const double shapeparam = (aspheric * aspheric + 0.75 * acylind * acylind) / (gyr_rad * gyr_rad);

  return shapeparam;
}

double Observer::CalcThickness(const dpdsystem& sDPD, const Parameter& param) {
  double up_pos = 0.0; int up_elem = 0;
  double dw_pos = 0.0; int dw_elem = 0;
  for (int i = 0; i < param.ampN; i++) {
    const bool chm  = sDPD.GetLipidChemConf(i);
    if (chm) {
      const int h_idx = sDPD.GetLipidElemIdx(Parameter::ALL_UNIT_N * i          );
      const int t_idx = sDPD.GetLipidElemIdx(Parameter::ALL_UNIT_N * (i + 1) - 1);
      
      const double3 h_pos = sDPD.pr[h_idx];
      const double3 t_pos = sDPD.pr[t_idx];
      double dry = h_pos.y - t_pos.y;
      dry -= param.L.y * std::round(dry * param.iL.y);
      if (dry > 0.0) {
	up_pos += h_pos.y;
	up_elem++;
      } else {
	dw_pos += h_pos.y;
	dw_elem++;
      }
    }
  }
  up_pos /= up_elem;
  dw_pos /= dw_elem;
  
  return std::fabs(up_pos - dw_pos);
}

double Observer::CalcOrientOrder(const dpdsystem& sDPD,const Parameter& param) {
  double orderS = 0.0;
  for (int i = 0; i < param.ampN; i++) {
    const bool chm = sDPD.GetLipidChemConf(i);
    if (chm) {
      const int h_idx = sDPD.GetLipidElemIdx(Parameter::ALL_UNIT_N * i          );
      const int t_idx = sDPD.GetLipidElemIdx(Parameter::ALL_UNIT_N * (i + 1) - 1);
      double3 dr = sDPD.pr[h_idx] - sDPD.pr[t_idx];
      dr.x       -= param.L.x * std::round(dr.x * param.iL.x);
      dr.y       -= param.L.y * std::round(dr.y * param.iL.y);
      dr.z       -= param.L.z * std::round(dr.z * param.iL.z);
      const double norm = std::sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
      dr /= norm;
      orderS    += dr.y * dr.y;
    }
  }
  orderS *= 1.5 / param.ampN;
  orderS -= 0.5;
  return orderS;
}

void Observer::CalcMembraneHeight(const dpdsystem& sDPD) {
  hei.assign(hei.size(), 0.0);
  hei_elem.assign(hei_elem.size(), 0);
  for (int i = 0; i < Parameter::sys_size; i++) {
    if (sDPD.prop[i] != 0) {
      const int hash = GenHash2d(sDPD.pr[i]);
      hei[hash] += sDPD.pr[i].y;
      hei_elem[hash]++;
    }
  }

  for (size_t i = 0; i < hei.size(); i++) {
    assert(hei_elem[i] != 0);
    hei[i] /= hei_elem[i];
  }
}

void Observer::DumpMacroVal(const dpdsystem& sDPD,const Parameter& param){
  const double kT = CalcKinTempera(sDPD);
  const double Dif = CalcDiffs(sDPD);
  const double thick = CalcThickness(sDPD, param);
  const double OrderS = CalcOrientOrder(sDPD, param);
  fprintf(fp[MACRO_VAL], "%.10g %.10g %.10g %.10g \n", kT, Dif, thick, OrderS);
}

void Observer::DumpPressure(const dpdsystem& sDPD, const Parameter& param, const double3& vil){
  const double3 kin_strs = std::accumulate(sDPD.pv, sDPD.pv + Parameter::sys_size, double3(0.0, 0.0, 0.0),
					   [](const double3& sum, const double3& val) {return sum + val;});
  const double3 P  = (vil + kin_strs) * param.iL.x * param.iL.y * param.iL.z ;
  const double Sig = (P.y - (P.x + P.z) * 0.5) * param.L.y;
  const double Gam = ( (P.x + P.y) * 0.5 - P.z ) * 0.5 * param.L.x * param.L.y;
  fprintf(fp[PRESSURE], "%.10g %.10g %.10g %.10g %.10g \n", P.x, P.y, P.z, Sig, Gam);
}

void Observer::DumpConfigTempera(const double configT) {
  fprintf(fp[CONFIG_TEMP], "%.10g \n", configT);
}

void Observer::DumpPrtclConfig(const dpdsystem &sDPD, const ChemInfo& cheminfo, const int time, FILE* fp) {
  static constexpr char atom_type[] = {
    'O', 'N', 'C', 'S'
  };

  fprintf(fp, "%d\n", Parameter::sys_size);
  fprintf(fp, "time %d\n", time);
  for (int i = 0; i < Parameter::sys_size; i++) {
    const double3 pos   = sDPD.pr[i];
    const double3 vel   = sDPD.pv[i];
    const par_prop prp  = sDPD.prop[i];
    const bool   p_chem = cheminfo.prtcl_chem[i];
    const int  part_idx = cheminfo.part_idx[i];
    const int  l_unit   = cheminfo.lipid_unit[i];
    const int  l_idx    = cheminfo.lipid_idx[i];

    //NOTE:
    //pos.x pos.y pos.z vel.x vel.y vel.z prop prtcl_chem
    //lipid_chem lipid_unit lipid_idx
    //part_idx
    if (prp == Water) {
      fprintf(fp, "%c %.15g %.15g %.15g %.15g %.15g %.15g %d %d %d %d %d %d \n",
	      atom_type[prp], pos.x, pos.y, pos.z, vel.x, vel.y, vel.z, prp, p_chem,
	      false, l_unit, l_idx,
	      part_idx);
    } else {
      bool chem = false;
      if (l_idx != -1) chem = sDPD.GetLipidChemConf(l_idx);
      fprintf(fp, "%c %.15g %.15g %.15g %.15g %.15g %.15g %d %d %d %d %d %d \n",
	      atom_type[prp], pos.x, pos.y, pos.z, vel.x, vel.y, vel.z, prp, p_chem,
	      chem, l_unit, l_idx,
	      part_idx);
    }
  }
}

void Observer::DumpLocalVal(const dpdsystem &sDPD, const Parameter& param) {
  const double3* t_r = sDPD.pr;
  const double3* t_v = sDPD.pv;

  loc_tempera.assign(loc_tempera.size(),0.0);
  loc_dense.assign(loc_dense.size(),0.0);
  loc_vel.assign(loc_vel.size(), double3(0.0));

  for (int i = 0; i < Parameter::sys_size; i++) {
    const double3 pos = t_r[i];
    const double3 vel = t_v[i];
    int hash = static_cast<int>(pos.z * param.i_grid_leng.z);
    if (hash == param.grid_numb[1]) hash--;
    loc_tempera[hash] += vel*vel;
    loc_dense[hash]   += 1.0;
    loc_vel[hash] += vel;
  }

  for (size_t i = 0; i < loc_tempera.size(); i++) {
    loc_vel[i] /= loc_dense[i];
    loc_tempera[i] /= loc_dense[i] * 3.0;
    loc_dense[i] /= param.grid_leng.z * param.L.x * param.L.y;
  }
  
  for (size_t i = 0; i < loc_tempera.size(); i++) {
    const double cur_z = (i + 0.5) * param.grid_leng.z;
    fprintf(fp[LOCAL_VAL], "%f %.15g %.15g %.15g %.15g %.15g\n", cur_z, loc_tempera[i], loc_dense[i], loc_vel[i].x, loc_vel[i].y, loc_vel[i].z);
  }
}

void Observer::DumpTranject(const dpdsystem& sDPD, const ChemInfo& cheminfo, const int time) {
  DumpPrtclConfig(sDPD, cheminfo, time, fp[PTCL_CONFIG]);
}

void Observer::DumpFinalConfig(const dpdsystem &sDPD, const ChemInfo& cheminfo, const int time) {
  DumpPrtclConfig(sDPD, cheminfo, time, fp[FINAL_CONFIG]);
}

void Observer::DumpMembHeight(const dpdsystem& sDPD, const Parameter& param, const int time) {
  CalcMembraneHeight(sDPD);
  fprintf(fp[HEIGHT_DIST] , "%f ", time * param.dt);
  for(size_t i = 0; i < hei.size(); i++) fprintf(fp[HEIGHT_DIST], "%.15g ", hei[i]);
  fprintf(fp[HEIGHT_DIST], "\n");
}
