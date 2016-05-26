#pragma once

#include <vector>
#include <array>
#include "mvector3.hpp"

//modes
//for debug
// #define NO_THERMO_STAT

//chem
// #define CHEM_MODE
#ifdef CHEM_MODE
#define INVERSE_CHEM
#endif

//boundary
// #define Y_REFLECT_BOUND

//bind potential
// #define ADD_BIND_POTENT

//position restraint
// #define ADD_POSRES

//
// #define CALC_HEIGHT

//visc
// #define VISC_KERN_SQRT

#define CALC_LOC_STRESS

#ifdef CALC_LOC_STRESS
// #define RELAXED_BASE_POS
#define CENTRAL_FORCE
#endif

enum par_prop {
  Water = 0, Hyphil, Hyphob,
  
  Numprop,
};

// Local stress type
enum {
  INTER_MOL = 0,
  BOND,
  ANGLE,
  // DIHEDRAL,
  KINETIC,
  
  NUM_LS_TYPE,
};

struct Interactions {
  double                cf_sigma[Numprop][Numprop];
  double                cf_gamma[Numprop][Numprop];
  double                cf_spring;
  double                cf_bend;
  double                cf_repul[Numprop][Numprop];
};

struct BindInfo {
  double3 bind_center;
  double  bind_radius;
  double  bind_coef;
};

struct CellList {
  int  *buck_elem, *buck_addrs, *next_dest;
  std::vector<int>* near_prtcl_idx;
  int *n_near_prtcl;
  int *prtcl_idx;
  int *prtcl_in_cell;
};

struct float_int{
  float dist;
  int idx;
};

struct ChemInfo{
  bool *prtcl_chem, *lipid_chem;
  int  *lipid_idx,  *part_idx;
  int  *lipid_unit;
  int  *lip_elem_idx;
  int  *head_elem_idx, *tail_elem_idx;
  float *near_water;
  float_int *near_info;
};

typedef std::array<std::array<double, 3>, 3> tensor3d;

inline const tensor3d operator + (const tensor3d& lhs, const tensor3d& rhs) {
  tensor3d ret = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
      ret[i][j] = lhs[i][j] + rhs[i][j];
    }
  return ret;
}

inline const tensor3d operator - (const tensor3d& lhs, const tensor3d& rhs) {
  tensor3d ret = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
      ret[i][j] = lhs[i][j] - rhs[i][j];
    }
  return ret;
}


#define CHECK_FILE_OPEN(fin)					\
  do {								\
    if (!fin) {							\
      std::cerr << "Cannot open " << fname << "." << std::endl;	\
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;	\
      std::exit(1);						\
    }								\
  } while (false)

#define CHECK_FILESTREAM_IS_OK(fin)				\
  do {								\
    if (fin.fail()) {						\
      std::cerr << "There is a problem in file stream.\n";	\
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;	\
      std::exit(1);						\
    }								\
  } while (false)

#define CHECK_FILE_IS_EOF(fin)					\
  do {								\
    double buf = 0.0;						\
    fin >> buf;							\
    if (!fin.eof()) {						\
      std::cerr << "is not eof\n.";				\
      std::cerr << __FILE__  << " " << __LINE__ << std::endl;	\
      std::exit(1);						\
    }								\
  } while (false)

#define CHECK_EQUATION(eqa, val)					\
  do {									\
    if (!(eqa)) {							\
      std::cerr << #eqa << " is not satisfied.\n";			\
      std::cerr << #val << " = " << val << std::endl;			\
      std::cerr << __FILE__ << " " << __LINE__ << std::endl;		\
      std::exit(1);							\
    }									\
  } while (false)

inline FILE* xfopen(const char* __restrict filename,
		    const char* __restrict  mode) {
  FILE* f = fopen(filename, mode);
  if (f == nullptr) {
    std::cerr << filename << ": Cannot open file\n.";
    std::exit(1);
  }
  return f;
}
