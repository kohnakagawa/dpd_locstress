#include "bucket_sorter.hpp"

B_sorter::B_sorter(const Parameter& param) {
  interact_list   = new std::array<int, 13> [param.all_grid];
  buf_double3	  = new double3  [Parameter::SYS_SIZE];
  buf_int	  = new int [Parameter::SYS_SIZE];
  buf_prop	  = new par_prop [Parameter::SYS_SIZE];
  buf_bool	  = new bool [Parameter::SYS_SIZE];
    
  GenInteractList(param);
}

B_sorter::~B_sorter() {
  delete [] interact_list;
  delete [] buf_double3;
  delete [] buf_int;
  delete [] buf_prop;
  delete [] buf_bool;
}

void B_sorter::GenInteractList(const Parameter& param) {
  int tar_grid = 0;
  for (int iz = 0; iz < param.grid_numb[2]; iz++) {
    for (int iy = 0; iy < param.grid_numb[1]; iy++) {
      for (int ix = 0; ix < param.grid_numb[0]; ix++) {
	int inr_grid = 0;
	for (int jz = -1; jz < 2; jz++) {
	  for (int jy = -1; jy < 2; jy++) {
	    for (int jx = -1; jx < 2; jx++) {
	      int q[] = {ix + jx, iy + jy, iz + jz};
	      ApplyPeriodicBoundary(q, param);
	      interact_list[tar_grid][inr_grid] = GenHash(q,param);
	      inr_grid++;
	      if (inr_grid == 13) goto OUT;
	    }
	  }
	}
      OUT:
	  
	tar_grid++;
      }
    }
  }
}

void B_sorter::BucketSort(dpdsystem& sDPD,
			  ChemInfo& cheminfo,
			  const Parameter& param,
			  const CellList& celllist) {
  CopyGather(sDPD.pr, buf_double3, Parameter::SYS_SIZE, celllist);
  CopyGather(sDPD.pv, buf_double3, Parameter::SYS_SIZE, celllist);
  CopyGather(sDPD.prop, buf_prop, Parameter::SYS_SIZE, celllist);
  CopyGather(sDPD.pv_bef, buf_double3, Parameter::SYS_SIZE, celllist);
  CopyGather(sDPD.delta_sumr, buf_double3, Parameter::SYS_SIZE, celllist);
  CopyGather(sDPD.force, buf_double3, Parameter::SYS_SIZE, celllist);
  CopyGather(sDPD.force_bef, buf_double3, Parameter::SYS_SIZE, celllist);
  CopyGather(cheminfo.prtcl_chem, buf_bool, Parameter::SYS_SIZE, celllist);
  CopyGather(cheminfo.lipid_idx, buf_int, Parameter::SYS_SIZE, celllist);
  CopyGather(cheminfo.part_idx, buf_int, Parameter::SYS_SIZE, celllist);
  CopyGather(cheminfo.lipid_unit, buf_int, Parameter::SYS_SIZE, celllist);

#ifdef ADD_POSRES
  CopyGather(sDPD.rest_on, buf_bool, Parameter::SYS_SIZE, celllist);
  CopyGather(sDPD.pr_base, buf_double3, Parameter::SYS_SIZE, celllist);
#endif
  
  for (int i = 0; i < Parameter::SYS_SIZE; i++)
    celllist.prtcl_idx[i] = i;

  for (int tar_grid = 0; tar_grid < param.all_grid; tar_grid++) {
    const int begid = celllist.buck_addrs[tar_grid    ];
    const int endid = celllist.buck_addrs[tar_grid + 1];
    for (int i = begid; i < endid; i++) {
      celllist.prtcl_in_cell[i] = tar_grid;
    }
  }
}

void B_sorter::MkPrtclIdx(const double3* pr,
			  const Parameter& param,
			  CellList& celllist) {
    //init bucket
    for (int i = 0; i < param.all_grid; i++) celllist.buck_elem[i] = 0;
    celllist.buck_addrs[0] = 0;
  
    //counting elem
    for (int i = 0; i < Parameter::SYS_SIZE; i++) {
      const int hash = GenHash(pr[i], param);
      celllist.buck_elem[hash]++;
      celllist.prtcl_in_cell[i] = hash;
    }

    //calc buck address
    for (int i = 0; i < param.all_grid; i++) 
      celllist.buck_addrs[i + 1] = celllist.buck_addrs[i] + celllist.buck_elem[i];

#ifdef DEBUG
    assert(celllist.buck_addrs[param.all_grid] == Parameter::SYS_SIZE);
#endif
  
    //calc destination
    for (int i = 0; i < Parameter::SYS_SIZE; i++) {
      const int	hash		   = celllist.prtcl_in_cell[i];
      const int	temp_idx	   = celllist.buck_addrs[hash];
      celllist.next_dest[i]	   = temp_idx;
      celllist.prtcl_idx[temp_idx] = i;
      celllist.buck_addrs[hash]++;
    }
  
    //recalc
    for (int i = 0; i < param.all_grid; i++) 
      celllist.buck_addrs[i] = celllist.buck_addrs[i] - celllist.buck_elem[i];
}

void B_sorter::ClearPrtclIdx(const Parameter& param, CellList& celllist) {
  for (int i = 0; i < param.all_grid; i++) {
    celllist.near_prtcl_idx[i].clear();
  }
}

void B_sorter::SetPrtclCellIdx(CellList& celllist, const Parameter& param) {
#pragma omp parallel for
  for (int tar_grid = 0; tar_grid < param.all_grid; tar_grid++) {
    const int begi = celllist.buck_addrs[tar_grid    ];
    const int endi = celllist.buck_addrs[tar_grid + 1];
    celllist.near_prtcl_idx[tar_grid].insert(celllist.near_prtcl_idx[tar_grid].end(),&celllist.prtcl_idx[begi],&celllist.prtcl_idx[endi]);
    for (int i = 0; i < 13; i++) {
      const int intr_grid = interact_list[tar_grid][i];
      const int begj	  = celllist.buck_addrs[intr_grid    ];
      const int endj	  = celllist.buck_addrs[intr_grid + 1];
      celllist.near_prtcl_idx[tar_grid].insert(celllist.near_prtcl_idx[tar_grid].end(), &celllist.prtcl_idx[begj], &celllist.prtcl_idx[endj]);
    }
  }
}

void B_sorter::SetPrtclCellIdxNoSTL(CellList& celllist, const Parameter& param) {
#pragma omp parallel for  
  for (int tar_grid = 0; tar_grid < param.all_grid; tar_grid++) {
    int n_intr_cell = 0;
    const int begi = celllist.buck_addrs[tar_grid    ];
    const int endi = celllist.buck_addrs[tar_grid + 1];
    for (int id = begi; id < endi; id++) {
      celllist.near_prtcl_idx[tar_grid][n_intr_cell] = celllist.prtcl_idx[id];
      n_intr_cell++;
    }
    for (int i = 0; i < 13; i++) {
      const int intr_grid = interact_list[tar_grid][i];
      const int begj	  = celllist.buck_addrs[intr_grid    ];
      const int endj	  = celllist.buck_addrs[intr_grid + 1];
      for (int id = begj; id < endj; id++) {
	celllist.near_prtcl_idx[tar_grid][n_intr_cell] = celllist.prtcl_idx[id];
	n_intr_cell++;
      }
    }
    celllist.n_near_prtcl[tar_grid] = n_intr_cell;
    
#ifdef DEBUG
    assert(n_intr_cell < Parameter::BUF_SIZE);
#endif    
  }
}

void B_sorter::CheckSorted(const dpdsystem& sDPD,
			   const Parameter& param,
			   const CellList& celllist) const {
  /*for (int i = 0; i < param.all_grid; i++) {
    for (int j = 0; j < 13; j++) {
    std::cout << interact_list[i][j] << std::endl;
    }
    }*/
  for (int tar_grid = 0; tar_grid < param.all_grid; tar_grid++) {
    for (int pi = celllist.buck_addrs[tar_grid]; pi < celllist.buck_addrs[tar_grid + 1]; pi++) {
      assert(pi < Parameter::SYS_SIZE);
      assert(celllist.prtcl_idx[pi] < Parameter::SYS_SIZE);
      const int hash = GenHash(sDPD.pr[celllist.prtcl_idx[pi]], param);
      if (hash != tar_grid) {
	std::cout << pi  << " " << hash << " " << tar_grid << std::endl;
	std::cout << sDPD.pr[pi] << std::endl;
      }
      assert(hash == tar_grid);
    }
  }
}

void B_sorter::CheckNearList(const dpdsystem& sDPD, const Parameter& param) const {
  std::ofstream fout("pair.txt");
  for (int i = 0; i < Parameter::SYS_SIZE - 1; i++) {
    for (int j = i + 1; j < Parameter::SYS_SIZE; j++) {
      double3 dr = sDPD.pr[i] - sDPD.pr[j];
      dr.x -= param.L.x * static_cast<int>(dr.x * param.ihL.x);
      dr.y -= param.L.y * static_cast<int>(dr.y * param.ihL.y);
      dr.z -= param.L.z * static_cast<int>(dr.z * param.ihL.z);

      const double dr2 = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
      if (dr2 < 1.0) {
	fout << i << " " << j << std::endl; 
      }
    }
  }
}

