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

class ConfigMaker {
  const char* cur_dir_ = nullptr;
  std::vector<PtclBuffer> ptcl_buffer_;

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
    cm_vel /= Parameter::SYS_SIZE;
    for (auto& ptcl : ptcl_buffer_) ptcl.vel -= cm_vel;
  }

  void SortPtclBuffer(const Parameter& param) {
    for (int i = 0; i < Parameter::SYS_SIZE; i++) {
      int key[3] = {-1};
      for (int j = 0; j < 3; j++) {
	key[j] = static_cast<int>(ptcl_buffer_[i].pos[j] * param.i_grid_leng[j]);
	if (key[j] == param.grid_numb[j]) key[j]--;
      }
      ptcl_buffer_[i].hash = key[0] + param.grid_numb[0] * (key[1] + key[2] * param.grid_numb[1]);
    }
    std::sort(ptcl_buffer_.begin(), ptcl_buffer_.end(), [](const PtclBuffer& i, const PtclBuffer& j) {return i.hash < j.hash;});
  }

  //NOTE: r = a base position n = the membrane normal vector
  //head case offset = 0                 end = param.head_part_n
  //tail case offset = param.head_part_n end = param.tail_part_n
  template<int offset, int end>
  void SetAmphilPartPos(const double3& r,
			const double3& n,
			const Parameter& param,
			int& phil_idx,
			int& phob_idx) {
    for (int unit = 0; unit < end; unit++) {
      double3 temp_pos = r + (unit + offset) * param.b_leng * n;
      ApplyPBC(temp_pos, param);
      if (unit + offset < Parameter::HYPHIL_N) {
	// hydrophilic part
	ptcl_buffer_[phil_idx].pos = temp_pos;
	ptcl_buffer_[phil_idx].prop = Hyphil;
	ptcl_buffer_[phil_idx].l_unit = unit + offset;
	phil_idx++;
      } else {
	//hydrophobic part
	ptcl_buffer_[phob_idx].pos = temp_pos;
	ptcl_buffer_[phob_idx].prop = Hyphob;
	ptcl_buffer_[phob_idx].l_unit = unit + offset;
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
      SetAmphilPartPos<0, Parameter::REAC_PART>(base, nv, param, phil_idx, phob_idx);
      SetAmphilPartPos<Parameter::REAC_PART, Parameter::TAIL_PART>(base, nv, param, phil_idx, phob_idx);
      added_num += Parameter::ALL_UNIT_N;
    }
  }

  void MakeSphCap(int& added_num, const Parameter& param) {
    // not implemented yet.
  }

  template<class prng_gen>
  void MakeRandom(prng_gen& prng,
		  std::uniform_real_distribution<>& uniform,
		  int& added_num,
		  const Parameter& param) {
    int phil_idx = 0, phob_idx = param.hN;
    for (int i = 0; i < param.ampN; i++) {
      const auto base	= GenRandVec(prng, uniform, param.L);
      auto nv		= Genrandvec(prng, uniform, double3(1.0, 1.0, 1.0));
      nv		= nv / std::sqrt(nv * nv);
      SetAmphilPartPos<0, Parameter::REAC_PART>(base, nv, param, phil_idx, phob_idx);
      SetAmphilPartPos<Parameter::REAC_PART, Parameter::TAIL_PART>(base, nv, param, phil_idx, phob_idx);
      added_num += Parameter::ALL_UNIT_N;
    }
  }

  // NOTE: valid for flat membrane
  template<class prng_gen>  
  void FillWaterPtclsRect(prng_gen& prng,
			  std::uniform_real_distribution<>& uniform,
			  int& added_num,
			  const Parameter& param) {
    const double3 org = std::accumulate(ptcl_buffer_.cbegin(), ptcl_buffer_.cbegin() + added_num, param.L,
					[](const PtclBuffer& sum, const PtclBuffer& val) {return (sum.pos < val.pos) ? sum : val;});
    const double3 top = std::accumulate(ptcl_buffer_.cbegin(), ptcl_buffer_.cbegin() + added_num, double3(0.0, 0.0, 0.0),
					[](const PtclBuffer& sum, const PtclBuffer& val) {return (sum.pos > val.pos) ? sum : val;});
    
    int water_idx = param.hN + param.bN;
    while (water_idx < Parameter::SYS_SIZE) {
      const auto base = Genrandvec(prng, uniform, param.L);
      if (base > org && base < top) {
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
      MakeFlatMembrane(param.L, 1, param.ampN, added_num, prng, uniform, param);
      FillWaterPtclsRect(prng, uniform, added_num, param);
    } else if (mode == "Random") {
      MakeRandom(prng, uniform, added_num, param);
      FillWaterPtclsRandom(prng, uniform, added_num, param);
    // } else if (mode == "SphCap") {
    } else {
      std::cerr << "" << std::endl;
      std::exit(1);
    }
  }
  
  void ApplyPBC(double3& r, const Parameter& param) {
    r.x -= std::floor(r.x * param.iL.x) * param.L.x;
    r.y -= std::floor(r.y * param.iL.y) * param.L.y;
    r.z -= std::floor(r.z * param.iL.z) * param.L.z;
  }
  
  void CheckConfigurationIsValid(const Parameter& param) {
    // range check
    for (const auto& ptcl : ptcl_buffer_) {
      for (int j = 0; j < 3; j++) {
	CHECK_EQUATION(ptcl.pos[j] <= param.L[j], ptcl.pos[j]);
	CHECK_EQUATION(ptcl.pos[j] >= 0.0, ptcl.pos[j]);
      }
    }
    
    // temperature check
    double3 vel2_sum = std::accumulate(ptcl_buffer_.cbegin(), ptcl_buffer_.cend(), double3(0.0, 0.0, 0.0),
				      [](const double3& sum, const PtclBuffer& val) {
					 return double3(sum.x + val.vel.x * val.vel.x,
							sum.y + val.vel.y * val.vel.y,
							sum.z + val.vel.y * val.vel.z);
				      });
    vel2_sum /= Parameter::SYS_SIZE;
    std::cerr << "Temperature = " << vel2_sum << std::endl;
  }
  
public:
  ConfigMaker(const char* cur_dir) {
    cur_dir_ = cur_dir;
  }

  void GenParticles(const Parameter& param) {
    const auto fname = std::string(cur_dir_) + "/init_config.xyz";
    std::ofstream fout(fname.c_str());
    
    // allocate ptcl buffer
    ptcl_buffer_.resize(Parameter::SYS_SIZE);

    // prng
    std::mt19937 mt_rnd(time(nullptr));
    std::normal_distribution<> normal(0.0, std::sqrt(param.tempera));
    std::cout << "Now create particle configuration.\n";
    std::cout << "Info:\n";
    std::cout << "# of particles is " << Parameter::SYS_SIZE << std::endl;
    
    std::cout << "Thermal velocity generation.\n";
    GenThermalVeloc(mt_rnd, normal);
    std::cout << "Done.\n";
    
    std::cout << "Configuration generation.\n";
    GenConfiguration(mt_rnd, param);
    std::cout << "Done.\n";
      
    std::cout << "In the end, we check the validity of this particle configuration.\n";
    CheckConfigurationIsValid(param);
    std::cout << "Done.\n";
  }
};

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage:\n";
    std::cerr << "argv[1] = target directory name.\n";
    std::exit(1);
  }
  
  Parameter param(argv[1]);
  param.LoadParam();
  param.LoadCheck();
  
  ConfigMaker cmaker(argv[1]);
  cmaker.GenParticles(param);
}
