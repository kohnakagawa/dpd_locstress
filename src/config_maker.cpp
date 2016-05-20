#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <ctime>
#include <numeric>
#include <algorithm>

#include "ptcl_buffer.hpp"
#include "parameter.hpp"

constexpr char PtclBuffer::atom_type[21];
int Parameter::sys_size;

class ConfigMaker {
  const char* cur_dir_ = nullptr;
  std::vector<PtclBuffer> ptcl_buffer_;
  int h_array_id = 0, t_array_id = 0; 

  template<class prng_gen>
  double3 GenRandVec(prng_gen& prng,
		     std::uniform_real_distribution<>& uniform,
		     const double3& val) {
    return double3(uniform(prng) * val.x, uniform(prng) * val.y, uniform(prng) * val.z);
  }
  
  template<class prng_gen>
  void GenThermalVeloc(prng_gen& prng,
		       std::normal_distribution<>& normal) {
    for (auto& ptcl : ptcl_buffer_) {
      for (int i = 0; i < 3; i++)
	ptcl.vel[i] = normal(prng);
    }

    // remove cm drift
    double3 cm_vel = std::accumulate(ptcl_buffer_.cbegin(), ptcl_buffer_.cend(), double3(0.0, 0.0, 0.0), 
				     [](const double3& sum, const PtclBuffer& val) {return sum + val.vel;});
    cm_vel /= Parameter::sys_size;
    for (auto& ptcl : ptcl_buffer_) ptcl.vel -= cm_vel;
  }

  void SortPtclBuffer(const Parameter& param) {
    for (int i = 0; i < Parameter::sys_size; i++) {
      int key[3] = {-1};
      for (int j = 0; j < 3; j++) {
	key[j] = static_cast<int>(ptcl_buffer_[i].pos[j] * param.i_grid_leng[j]);
	if (key[j] == param.grid_numb[j]) key[j]--;
      }
      ptcl_buffer_[i].hash = key[0] + param.grid_numb[0] * (key[1] + key[2] * param.grid_numb[1]);
    }
    std::sort(ptcl_buffer_.begin(), ptcl_buffer_.end(), [](const PtclBuffer& i, const PtclBuffer& j) {return i.hash < j.hash;});
  }

  void SetPartialLipidId(const int lipid_id,
			 int& base_idx,
			 const int unit_leng,
			 const int num_init_bond) {
    for (int unit = 0; unit < unit_leng; unit++) {
      ptcl_buffer_[base_idx].l_idx = lipid_id;
      if (lipid_id < num_init_bond) {
	ptcl_buffer_[base_idx].l_chem = true;
	ptcl_buffer_[base_idx].p_chem = true;
      }
      base_idx++;
    }
  }

  void SetLipidId(const Parameter& param,
		  const int num_init_bond) {
    CHECK_EQUATION(num_init_bond <= param.ampN, num_init_bond);
    // hydrophilic: base_idx[0] hydrophobic: base_idx[1]
    int base_idx[2] = {0, param.hN};
  
    for (int lid = 0; lid < num_init_bond; lid++) {
      SetPartialLipidId(lid, base_idx[0], Parameter::HYPHIL_N, num_init_bond);
      SetPartialLipidId(lid, base_idx[1], Parameter::HYPHOB_N, num_init_bond);
      // setlipidchem(lid, true);
    }

    if (Parameter::HYPHIL_N >= Parameter::REAC_PART) {
      if (param.tailN <= param.headN) base_idx[0] += Parameter::REAC_PART * (param.headN - num_init_bond);
    } else {
      if (param.tailN <= param.headN) base_idx[1] += (Parameter::REAC_PART - Parameter::HYPHIL_N) * (param.headN - num_init_bond);
    }
  
    for (int lid = num_init_bond; lid < param.ampN; lid++) {
      if (Parameter::HYPHIL_N >= Parameter::REAC_PART) {
	if (param.tailN <= param.headN) {
	  // tail
	  SetPartialLipidId(lid, base_idx[0], Parameter::HYPHIL_N - Parameter::REAC_PART, num_init_bond);
	  SetPartialLipidId(lid, base_idx[1], Parameter::HYPHOB_N, num_init_bond);
	} else {
	  // head
	  SetPartialLipidId(lid, base_idx[0], Parameter::REAC_PART, num_init_bond);
	}
      } else {
	if (param.tailN <= param.headN) {
	  // tail
	  SetPartialLipidId(lid, base_idx[1], Parameter::TAIL_PART, num_init_bond);
	} else {
	  // head
	  SetPartialLipidId(lid, base_idx[0], Parameter::HYPHIL_N, num_init_bond);
	  SetPartialLipidId(lid, base_idx[1], Parameter::REAC_PART - Parameter::HYPHIL_N, num_init_bond);
	}
      }
    }
  }

  //NOTE: r = a base position n = the membrane normal vector
  //head case offset = 0                 end = param.head_part_n
  //tail case offset = param.head_part_n end = param.tail_part_n
  template<int offset, int end>
  void SetAmphilPartPos(const double3& r,
			const double3& n,
			const Parameter& param,
			int& phil_idx,
			int& phob_idx,
			const bool is_head) {
    for (int unit = 0; unit < end; unit++) {
      double3 temp_pos = r + (unit + offset) * param.b_leng * n;
      ApplyPBC(temp_pos, param);
      if (unit + offset < Parameter::HYPHIL_N) {
	// hydrophilic part
	ptcl_buffer_[phil_idx].pos = temp_pos;
	ptcl_buffer_[phil_idx].prop = Hyphil;
	ptcl_buffer_[phil_idx].l_unit = unit + offset;

	if (is_head) {
	  ptcl_buffer_[phil_idx].p_idx = h_array_id / Parameter::REAC_PART;
	  h_array_id++;
	} else {
	  ptcl_buffer_[phil_idx].p_idx = t_array_id / Parameter::TAIL_PART;
	  t_array_id++;
	}
	phil_idx++;
      } else {
	//hydrophobic part
	ptcl_buffer_[phob_idx].pos = temp_pos;
	ptcl_buffer_[phob_idx].prop = Hyphob;
	ptcl_buffer_[phob_idx].l_unit = unit + offset;

	if (is_head) {
	  ptcl_buffer_[phob_idx].p_idx = h_array_id / Parameter::REAC_PART;
	  h_array_id++;
	} else {
	  ptcl_buffer_[phob_idx].p_idx = t_array_id / Parameter::TAIL_PART;
	  t_array_id++;
	}
	phob_idx++;
      }
    }
  }

  void SetWaterPos(const double3& r, int& water_idx) {
    ptcl_buffer_[water_idx].pos		= r;
    ptcl_buffer_[water_idx].prop	= Water;
    ptcl_buffer_[water_idx].l_unit	= -1;
    water_idx++;
  }
  
  template<class prng_gen>
  void MakeFlatMembrane(const double3& len,
			const int axis,
			const int amp_num,
			int& added_num,
			prng_gen& prng,
			std::uniform_real_distribution<>& uniform,
			const Parameter& param) {
    CHECK_EQUATION(axis >= 0, axis);
    CHECK_EQUATION(axis < 3, axis);
    
    int phil_idx = 0, phob_idx = param.hN;
    
    bool flap = true;
    const double eps = 1.0e-5;
    double3 nv(0.0, 0.0, 0.0);
    for (int i = 0; i < amp_num; i++) {
      double3 base = GenRandVec(prng, uniform, param.L);
      const double sign = flap ? 1.0 : -1.0;
      base[axis] = 0.5 * len[axis] + ((Parameter::ALL_UNIT_N - 1) * Parameter::b_leng + eps) * sign;
      nv[axis]   = -sign;
      SetAmphilPartPos<0, Parameter::REAC_PART>(base, nv, param, phil_idx, phob_idx, true);
      SetAmphilPartPos<Parameter::REAC_PART, Parameter::TAIL_PART>(base, nv, param, phil_idx, phob_idx, false);
      added_num += Parameter::ALL_UNIT_N;
      flap ^= true;
    }
  }

  // void MakeSphCap(int& added_num, const Parameter& param) {
  //   // not implemented yet.
  // }

  template<class prng_gen>
  void MakeRandom(prng_gen& prng,
		  std::uniform_real_distribution<>& uniform,
		  int& added_num,
		  const Parameter& param) {
    int phil_idx = 0, phob_idx = param.hN;
    for (int i = 0; i < param.ampN; i++) {
      const auto base	= GenRandVec(prng, uniform, param.L);
      auto nv		= GenRandVec(prng, uniform, double3(1.0, 1.0, 1.0));
      const auto norm   = std::sqrt(nv * nv);
      nv /= norm;
      SetAmphilPartPos<0, Parameter::REAC_PART>(base, nv, param, phil_idx, phob_idx, true);
      SetAmphilPartPos<Parameter::REAC_PART, Parameter::TAIL_PART>(base, nv, param, phil_idx, phob_idx, false);
      added_num += Parameter::ALL_UNIT_N;
    }
  }

  // NOTE: valid for flat membrane
  template<class prng_gen>  
  void FillWaterPtclsRect(prng_gen& prng,
			  std::uniform_real_distribution<>& uniform,
			  int& added_num,
			  const int axis,
			  const Parameter& param) {
    const double org = std::accumulate(ptcl_buffer_.cbegin(), ptcl_buffer_.cbegin() + added_num, param.L[axis],
				       [=](const double sum, const PtclBuffer& val) {return (sum < val.pos[axis]) ? sum : val.pos[axis];});
    const double top = std::accumulate(ptcl_buffer_.cbegin(), ptcl_buffer_.cbegin() + added_num, 0.0,
				       [=](const double sum, const PtclBuffer& val) {return (sum > val.pos[axis]) ? sum : val.pos[axis];});
    int water_idx = param.hN + param.bN;
    while (water_idx < Parameter::sys_size) {
      const auto base = GenRandVec(prng, uniform, param.L);
      if (!(base[axis] > org && base[axis] < top)) {
	SetWaterPos(base, water_idx);
	added_num++;
      }
    }
  }

  template<class prng_gen>
  void FillWaterPtclsRandom(prng_gen& prng,
			    std::uniform_real_distribution<>& uniform,
			    int& added_num,
			    const Parameter& param) {
    int water_idx = param.hN + param.bN;
    for (int i = 0; i < param.wN; i++) {
      SetWaterPos(GenRandVec(prng, uniform, param.L), water_idx);
      added_num++;
    }
  }

  template<class prng_gen>
  void GenConfiguration(prng_gen& prng, const Parameter& param) {
    std::cout << "Choose membrane configuration.\n";
    std::cout << "[Flat/Random]\n";
    std::string mode;
    std::cin >> mode;

    std::uniform_real_distribution<> uniform(0.0, 1.0);    
    int added_num = 0;
    if (mode == "Flat") {
      const int axis = 1;
      MakeFlatMembrane(param.L, axis, param.ampN, added_num, prng, uniform, param);
      FillWaterPtclsRect(prng, uniform, added_num, axis, param);
    } else if (mode == "Random") {
      MakeRandom(prng, uniform, added_num, param);
      FillWaterPtclsRandom(prng, uniform, added_num, param);
    // } else if (mode == "SphCap") {
    } else {
      std::cerr << "Unknown execution mode.\n";
      std::exit(1);
    }
  }
  
  void ApplyPBC(double3& r, const Parameter& param) {
    r.x -= std::floor(r.x * param.iL.x) * param.L.x;
    r.y -= std::floor(r.y * param.iL.y) * param.L.y;
    r.z -= std::floor(r.z * param.iL.z) * param.L.z;
  }
  
  void CheckConfigurationIsValid(const Parameter& param) {
    // position range check
    std::cout << "Range check.\n";
    for (const auto& ptcl : ptcl_buffer_) {
      for (int j = 0; j < 3; j++) {
	CHECK_EQUATION(ptcl.pos[j] <= param.L[j], ptcl.pos[j]);
	CHECK_EQUATION(ptcl.pos[j] >= 0.0, ptcl.pos[j]);
      }
    }
    std::cout << "Done.\n\n";
    
    // temperature check
    std::cout << "Temperature check.\n";
    double3 vel2_sum = std::accumulate(ptcl_buffer_.cbegin(), ptcl_buffer_.cend(), double3(0.0, 0.0, 0.0),
				       [](const double3& sum, const PtclBuffer& val) {
					 return double3(sum.x + val.vel.x * val.vel.x,
							sum.y + val.vel.y * val.vel.y,
							sum.z + val.vel.z * val.vel.z);
				       });
    vel2_sum /= Parameter::sys_size;
    
    std::cout << "Temperature = " << vel2_sum << std::endl;
    std::cout << "Reference temperature = " << param.tempera << std::endl;
    std::cout << "Done.\n\n";

    // topology check
    std::cout << "Topology check.\n";
    std::array<int, Parameter::ALL_UNIT_N> num_in_unit = {0};
    std::vector<int> num_in_pidx, num_in_lidx;
    num_in_pidx.resize(param.ampN, 0); num_in_lidx.resize(param.ampN, 0);
    
    for (const auto& ptcl : ptcl_buffer_) {
      if (ptcl.prop == Water) {
	CHECK_EQUATION(ptcl.p_chem == false, ptcl.p_chem);
	CHECK_EQUATION(ptcl.l_chem == false, ptcl.l_chem);
	CHECK_EQUATION(ptcl.l_unit == -1, ptcl.l_unit);
	CHECK_EQUATION(ptcl.l_idx  == -1, ptcl.l_idx);
	CHECK_EQUATION(ptcl.p_idx  == -1, ptcl.p_idx);
      } else {
	CHECK_EQUATION(ptcl.p_chem == true, ptcl.p_chem);
	CHECK_EQUATION(ptcl.l_chem == true, ptcl.l_chem);
	CHECK_EQUATION(ptcl.l_unit != -1, ptcl.l_unit);
	CHECK_EQUATION(ptcl.l_idx  != -1, ptcl.l_idx);
	CHECK_EQUATION(ptcl.p_idx  != -1, ptcl.p_idx);
	
	CHECK_EQUATION(ptcl.l_idx >= 0, ptcl.l_idx);
	CHECK_EQUATION(ptcl.l_idx < param.ampN, ptcl.l_idx);
	CHECK_EQUATION(ptcl.p_idx >= 0, ptcl.p_idx);
	CHECK_EQUATION(ptcl.p_idx < param.ampN, ptcl.p_idx);

	num_in_unit[ptcl.l_unit]++;
	num_in_pidx[ptcl.p_idx]++;
	num_in_lidx[ptcl.l_idx]++;
      }
    }
    
    CHECK_EQUATION(h_array_id == param.headN * Parameter::REAC_PART, h_array_id);
    CHECK_EQUATION(t_array_id == param.tailN * Parameter::TAIL_PART, t_array_id);
    
    for (int i = 0; i < Parameter::ALL_UNIT_N; i++)
      CHECK_EQUATION(num_in_unit[i] == param.ampN, num_in_unit[i]);

    for (int i = 0; i < param.ampN; i++) {
      CHECK_EQUATION(num_in_pidx[i] == Parameter::ALL_UNIT_N, num_in_pidx[i]);
      CHECK_EQUATION(num_in_lidx[i] == Parameter::ALL_UNIT_N, num_in_lidx[i]);
    }
    std::cout << "Done.\n\n";
  }

  void WritePtclData(std::ostream& ost) const {
    ost << ptcl_buffer_.size() << std::endl;
    ost << "time 0\n";
    for (const auto& ptcl : ptcl_buffer_) ptcl.WriteOstream(ost);
  }
  
public:
  explicit ConfigMaker(const char* cur_dir) {
    cur_dir_ = cur_dir;
  }

  void GenParticles(const Parameter& param) {
    const auto fname = std::string(cur_dir_) + "/init_config.xyz";
    std::ofstream fout(fname.c_str());
    
    // allocate ptcl buffer
    ptcl_buffer_.resize(Parameter::sys_size);

    // prng
    std::mt19937 mt_rnd(time(nullptr));
    std::normal_distribution<> normal(0.0, std::sqrt(param.tempera));
    std::cout << "Now create particle configuration.\n\n";
    
    std::cout << "Info:\n";
    param.DumpAllParam(std::cout);
    std::cout << std::endl;
    
    std::cout << "Thermal velocity generation.\n";
    GenThermalVeloc(mt_rnd, normal);
    std::cout << "Done.\n\n";
    
    std::cout << "Configuration generation.\n";
    GenConfiguration(mt_rnd, param);
    std::cout << "Done.\n\n";

    std::cout << "Set lipid id.\n";
    SetLipidId(param, param.ampN);
    std::cout << "Done.\n\n";
      
    std::cout << "In the end, we check the validity of this particle configuration.\n";
    CheckConfigurationIsValid(param);
    std::cout << "Done.\n\n";

    std::cout << "Will write configuration to " << fname << ".\n";
    WritePtclData(fout);
    std::cout << "Done.\n\n";
  }
};

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage:\n";
    std::cerr << "$" << argv[0] << " target directory .\n";
    std::cerr << "argv[1] = target directory name.\n";
    std::exit(1);
  }

  char* const dir_name = argv[1];
  
  Parameter param(dir_name);
  param.LoadParam();
  param.LoadCheck();
  
  ConfigMaker cmaker(dir_name);
  cmaker.GenParticles(param);

  std::cout << "Particle configuration is generaged at " << argv[1] << std::endl;
}
