#pragma once

#ifdef __INTEL_COMPILER
#include <mkl_lapacke.h>
template<int nRows, int nCols, int nRHS, int layout = LAPACK_COL_MAJOR>
void call_mkl_dgelsd(std::array<double, nRows * nCols>& D,
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
#endif

// F77 style
extern "C" {
  void dgelsd_(int* m, int* n, int* nrhs, double* a, int* lda,
	       double* b, int* ldb, double* s, double* rcond, int* rank,
	       double* work, int* lwork, int* iwork, int* info);
}

// Column-major order
template<int nRows, int nCols, int nRHS, int layout = 0>
void call_f77_dgelsd(std::array<double, nRows * nCols>& D,
		     std::array<double, nRows>& b,
		     double rcond = 1.0e-12) {
  std::array<double, nRows> s; s.fill(0.0);

  // From https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgelsd_ex.c.htm
  int nrows = nRows, ncols = nCols, nrhs = nRHS, info = -1, lwork = -1, rank = -1;
  int leading = (nrows > ncols) ? nrows : ncols;
  int iwork[3 * nRows * 0 + 11 * nCols];
  double wkopt = 0.0;

  dgelsd_(&nrows, &ncols, &nrhs, &D[0], &nrows,
	  &b[0], &leading, &s[0], &rcond, &rank,
	  &wkopt, &lwork, iwork, &info);
  lwork = static_cast<int>(wkopt);
  double* work = new double [lwork];
  dgelsd_(&nrows, &ncols, &nrhs, &D[0], &nrows,
	  &b[0], &leading, &s[0], &rcond, &rank,
	  work, &lwork, iwork, &info);
  if (info > 0) {
    std::cerr << "The algorithm computing SVD failed to converge\n";
    std::cerr << "The least squares solution could not be computed.\n";
    std::exit(1);
  }
  delete [] work;
}

