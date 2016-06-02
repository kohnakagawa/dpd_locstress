#pragma once

#ifdef __INTEL_COMPILER
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
template<int nRows, int nCols, int nRHS, int layout = LAPACK_COL_MAJOR>
void call_dgelsd(std::array<double, nRows * nCols>& D,
		 std::array<double, nRows>& b,
		 double rcond = 1.0e-12) {
  constexpr int leading = (nRows > nCols) ? nRows : nCols;

  std::array<double, nRows> s; s.fill(0.0);
  int rank = -1;
  
  const auto info = LAPACKE_dgelsd(layout,
				   nRows, nCols, nRHS, &D[0], nRows,
				   &b[0], leading, &s[0], rcond, &rank);
  if (info > 0) {
    std::cerr << "The algorithm computing SVD failed to converge\n";
    std::cerr << "The least squares solution could not be computed.\n";
    std::exit(1);
  }
}
