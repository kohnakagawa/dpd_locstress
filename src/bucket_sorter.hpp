#pragma once

#include "parameter.hpp"
#include "dpdsystem.hpp"
#include <array>

class B_sorter {
  std::array<int, 13>*  interact_list = nullptr;
  double3*              buf_double3   = nullptr;
  int*                  buf_int	      = nullptr;
  par_prop*             buf_prop      = nullptr;
  bool*                 buf_bool      = nullptr;

  template <typename T>
  void CopyGather(T* __restrict dat,
		  T* __restrict buf,
		  int elem,
		  const CellList& __restrict celllist) {
#pragma omp parallel for
    for (int i = 0; i < elem; i++)
      buf[i] = dat[i];
#ifdef DEBUG
    for (int i = 0; i < elem; i++)
      assert(celllist.next_dest[i] < Parameter::sys_size);
#endif
#pragma omp parallel for
    for (int i = 0; i < elem; i++)
      dat[celllist.next_dest[i]] = buf[i];
  }

  int GenHash(const int* q, const Parameter& param) const {
    const int hash = q[0] + param.grid_numb[0] * (q[1] + q[2] * param.grid_numb[1]);
#ifdef DEBUG
    assert(0 <= hash);
    assert(hash < param.all_grid);    
#endif
    return hash;
  }

  void ApplyPeriodicBoundary(int* q, const Parameter& param) const {
    if (q[0]<0 || q[0] >= param.grid_numb[0]) 
      q[0] -= static_cast<int>(q[0] * param.grid_numb[0] / std::abs(q[0]));

    if (q[1]<0 || q[1] >= param.grid_numb[1]) 
      q[1] -= static_cast<int>(q[1] * param.grid_numb[1] / std::abs(q[1]));

    if (q[2]<0 || q[2] >= param.grid_numb[2]) 
      q[2] -= static_cast<int>(q[2] * param.grid_numb[2] / std::abs(q[2]));
  }

  int GenHashForReflectBoundary(int* q, const Parameter& param, const int axis) const {
    if (q[axis] < 0 || q[axis] >= param.grid_numb[axis]) return -1;

    for (int i = 0; i < 3; i++) {
      if (i != axis) {
	if (q[i] < 0 || q[i] >= param.grid_numb[i])
	  q[i] -= static_cast<int>(q[i] * param.grid_numb[i] / std::abs(q[i]));
      }
    }
    
    return GenHash(q,param);
  }
  
  void GenInteractList(const Parameter& param);
public:
  enum {SORT_FREQ = 50};
  explicit B_sorter(const Parameter& param);
  ~B_sorter();
  void BucketSort(dpdsystem& sDPD, ChemInfo& cheminfo, const Parameter& param, const CellList& celllist);
  void MkPrtclIdx(const double3* pr, const Parameter& param, CellList& celllist);
  void ClearPrtclIdx(const Parameter& param, CellList& celllist);
  void SetPrtclCellIdx(CellList& celllist, const Parameter& param);
  void SetPrtclCellIdxNoSTL(CellList& celllist, const Parameter& param);

  void CheckSorted(const dpdsystem& sDPD, const Parameter& param, const CellList& celllist) const; 
  void CheckNearList(const dpdsystem& sDPD, const Parameter& param) const;

  static inline int GenHash(const double3 &r, const Parameter& param) {
    int i = static_cast<int>(r.x * param.i_grid_leng.x);
    int j = static_cast<int>(r.y * param.i_grid_leng.y);
    int k = static_cast<int>(r.z * param.i_grid_leng.z);

    if (i == param.grid_numb[0]) i--;
    if (j == param.grid_numb[1]) j--;
    if (k == param.grid_numb[2]) k--;
  
    const int hash = i + param.grid_numb[0] * (j + k * param.grid_numb[1]);
#ifdef DEBUG
    if(hash < 0 || hash > param.all_grid) std::cout << hash << std::endl;
    if(hash < 0 || hash > param.all_grid) std::cout << r << std::endl;
    assert(0 <= hash);
    assert(hash < param.all_grid);
#endif
    return hash;
  }
};
