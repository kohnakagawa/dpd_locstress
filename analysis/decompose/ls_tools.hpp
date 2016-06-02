#include "../../src/parameter.hpp"
#include <array>
#include <iomanip>

struct LSParam {
  std::array<double, 2> bleng = {{std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN()}};
  // NOTE: degree
  double bend_angle = std::numeric_limits<double>::signaling_NaN();
  double dip_angle  = std::numeric_limits<double>::signaling_NaN();
  
  void LoadParam(const char* cur_dir) {
    const std::string fname = std::string(cur_dir) + "/lsparam.txt";
    std::ifstream fin(fname.c_str());
    CHECK_FILE_OPEN(fin);
    fin >> bleng[0] >> bleng[1] >> bend_angle >> dip_angle;
    CHECK_FILESTREAM_IS_OK(fin);
    CHECK_FILE_IS_EOF(fin);
    
    // check the validity of this input parameter.
    CHECK_EQUATION(bend_angle > 0.0, bend_angle);
    CHECK_EQUATION(bend_angle < 180.0, bend_angle);
    CHECK_EQUATION(dip_angle > 0.0, dip_angle);
    CHECK_EQUATION(dip_angle < 180.0, dip_angle);
    CHECK_EQUATION(std::isfinite(bleng[0]), bleng[0]);
    CHECK_EQUATION(std::isfinite(bleng[1]), bleng[1]);

    bend_angle *= M_PI / 180.0;
    dip_angle  *= M_PI / 180.0;
  }
};

// LAPACK routine
extern "C" {
  void dgelsd_(int* m, int* n, int* nrhs, double* a, int* lda,
	       double* b, int* ldb, double* s, double* rcond, int* rank,
	       double* work, int* lwork, int* iwork, int* info);
}

template<int nRows, int nCols, int nRHS>
void call_dgelsd(double* D,
		 double* b,
		 double rcond = 1.0e-12) {
  std::array<double, nCols> s; s.fill(0.0);

  // From https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgelsd_ex.c.htm
  int nrows = nRows, ncols = nCols, nrhs = nRHS, info = -1, lwork = -1, rank = -1;
  int leading = (nrows > ncols) ? nrows : ncols;
  
  const int smlsiz = 25;
  const int nlvl = std::max(0, static_cast<int>(std::log2(std::min(nRows, nCols) / (smlsiz + 1))) + 1);
  const int iwork_size = 3 * std::min(nRows, nCols) * nlvl + 11 * std::min(nRows, nCols);
  int* iwork = new int [iwork_size];
  double wkopt = 0.0;

  dgelsd_(&nrows, &ncols, &nrhs, &D[0], &nrows,
	  &b[0], &leading, &s[0], &rcond, &rank,
	  &wkopt, &lwork, iwork, &info);
  lwork = static_cast<int>(wkopt);
  
  std::cout << "lwork " << lwork << std::endl;
  double* work = new double [lwork];
  dgelsd_(&nrows, &ncols, &nrhs, &(D[0]), &nrows,
	  &b[0], &leading, &s[0], &rcond, &rank,
	  work, &lwork, iwork, &info);
  if (info > 0) {
    std::cerr << "The algorithm computing SVD failed to converge\n";
    std::cerr << "The least squares solution could not be computed.\n";
    std::exit(1);
  }
  delete [] iwork;
  delete [] work;
}

double3 get_relaxed_pos(const double3& r1, const double3& r2,
			const double3& F1, const double3& F2) {
  constexpr int nRows = 3, nCols = 2, nRHS = 1;
  std::array<double, nRows * nCols> D;
  std::array<double, nRows> b;
  D[nRows * 0 + 0] = -F1.x; D[nRows * 1 + 0] = F2.x;
  D[nRows * 0 + 1] = -F1.y; D[nRows * 1 + 1] = F2.y;
  D[nRows * 0 + 2] = -F1.z; D[nRows * 1 + 2] = F2.z;
  b[0] = r1.x - r2.x; b[1] = r1.y - r2.y; b[2] = r1.z - r2.z;
    
  call_dgelsd<nRows, nCols, nRHS>(&D[0], &b[0]);
    
  const double3 cross_point = r1 + b[0] * F1;
    
  return cross_point;
}

void decompose_4n(const double3& r1,
		  const double3& r2,
		  const double3& r3,
		  const double3& r4,
		  const double3& r5,
		  const double3& F1,
		  const double3& F2,
		  const double3& F3,
		  const double3& F4,
		  double& f_decomp_err,
		  double& stress,
		  double3* dF) {
  constexpr int nRows = 15, nCols = 9, nRhs = 1;
  std::vector<double> D(nRows * nCols, 0.0);
  std::vector<double> b(nRows, 0.0);
  
  std::fill(D.begin(), D.end(), 0.0);
  std::fill(b.begin(), b.end(), 0.0);
  
  // do not consider PBC
  auto dr21 = r1 - r2;
  auto dr31 = r1 - r3;
  auto dr41 = r1 - r4;
  auto dr23 = r3 - r2;
  auto dr34 = r4 - r3;
  auto dr51 = r1 - r5;
  auto dr52 = r2 - r5;
  auto dr53 = r3 - r5;
  auto dr54 = r4 - r5;
  
  D[nRows * 0 + 0] = dr21.x; D[nRows * 0 + 3] = -dr21.x;
  D[nRows * 0 + 1] = dr21.y; D[nRows * 0 + 4] = -dr21.y;
  D[nRows * 0 + 2] = dr21.z; D[nRows * 0 + 5] = -dr21.y;
  
  D[nRows * 1 + 0] = dr31.x; D[nRows * 1 + 6] = -dr31.x;
  D[nRows * 1 + 1] = dr31.y; D[nRows * 1 + 7] = -dr31.y;
  D[nRows * 1 + 2] = dr31.z; D[nRows * 1 + 8] = -dr31.z;
  
  D[nRows * 2 + 0] = dr41.x; D[nRows * 2 + 9 ] = -dr41.x;
  D[nRows * 2 + 1] = dr41.y; D[nRows * 2 + 10] = -dr41.y;
  D[nRows * 2 + 2] = dr41.z; D[nRows * 2 + 11] = -dr41.z;
  
  D[nRows * 3 + 3] = -dr23.x; D[nRows * 3 + 6] = dr23.x;
  D[nRows * 3 + 4] = -dr23.y; D[nRows * 3 + 7] = dr23.y;
  D[nRows * 3 + 5] = -dr23.z; D[nRows * 3 + 8] = dr23.z;
  
  D[nRows * 4 + 6] = -dr34.x; D[nRows * 4 + 9]  = dr34.x;
  D[nRows * 4 + 7] = -dr34.y; D[nRows * 4 + 10] = dr34.y;
  D[nRows * 4 + 8] = -dr34.z; D[nRows * 4 + 11] = dr34.z;
  
  D[nRows * 5 + 0] = dr51.x; D[nRows * 5 + 12] = -dr51.x;
  D[nRows * 5 + 1] = dr51.y; D[nRows * 5 + 13] = -dr51.y;
  D[nRows * 5 + 2] = dr51.z; D[nRows * 5 + 14] = -dr51.z;

  D[nRows * 6 + 3] = dr52.x; D[nRows * 6 + 12] = -dr52.x;
  D[nRows * 6 + 4] = dr52.y; D[nRows * 6 + 13] = -dr52.y;
  D[nRows * 6 + 5] = dr53.z; D[nRows * 6 + 14] = -dr52.z;

  D[nRows * 7 + 6] = dr53.x; D[nRows * 7 + 12] = -dr53.x;
  D[nRows * 7 + 7] = dr53.y; D[nRows * 7 + 13] = -dr53.y;
  D[nRows * 7 + 8] = dr53.z; D[nRows * 7 + 14] = -dr53.z;

  D[nRows * 8 + 9 ] = dr54.x; D[nRows * 8 + 12] = -dr54.x;
  D[nRows * 8 + 10] = dr54.y; D[nRows * 8 + 13] = -dr54.y;
  D[nRows * 8 + 11] = dr54.z; D[nRows * 8 + 14] = -dr54.z;
  
  b[0] = F1.x; b[1] = F1.y; b[2] = F1.z;
  b[3] = F2.x; b[4] = F2.y; b[5] = F2.z;
  b[6] = F3.x; b[7] = F3.y; b[8] = F3.z;
  b[9] = F4.x; b[10] = F4.y; b[11] = F4.z;
  b[12] = 0.0; b[13] = 0.0; b[14] = 0.0;

  call_dgelsd<nRows, nCols, nRhs>(&D[0], &b[0]);

  const double3 dF21 = dr21 * b[0];
  const double3 dF31 = dr31 * b[1];
  const double3 dF41 = dr41 * b[2];
  const double3 dF23 = dr23 * b[3];
  const double3 dF34 = dr34 * b[4];
  const double3 dF51 = dr51 * b[5];
  const double3 dF52 = dr52 * b[6];
  const double3 dF53 = dr53 * b[7];
  const double3 dF54 = dr54 * b[8];

  // check err
  const auto F1_d = dF21 + dF31 + dF41 + dF51;
  const auto F2_d = dF52 - dF21 - dF23;
  const auto F3_d = dF23 + dF53 - dF31 - dF34;
  const auto F4_d = dF34 + dF54 - dF41;
  
  f_decomp_err += (F1_d - F1) * (F1_d - F1);
  f_decomp_err += (F2_d - F2) * (F2_d - F2);
  f_decomp_err += (F3_d - F3) * (F3_d - F3);
  f_decomp_err += (F4_d - F4) * (F4_d - F4);
  
  dF[0] = dF21;
  dF[1] = dF41;
  dF[2] = dF41;
  dF[3] = dF23;
  dF[4] = dF34;
  dF[5] = dF51;
  dF[6] = dF52;
  dF[7] = dF53;
  dF[8] = dF54;

  stress = 
    dF21.norm2() * dr21.norm2() +
    dF31.norm2() * dr31.norm2() +
    dF41.norm2() * dr41.norm2() +
    dF23.norm2() * dr23.norm2() +
    dF34.norm2() * dr34.norm2() +
    dF51.norm2() * dr51.norm2() +
    dF52.norm2() * dr52.norm2() +
    dF53.norm2() * dr53.norm2() +
    dF54.norm2() * dr54.norm2();
}

// from src/force_calculator.hpp
void decompose_3n(const double3& r1,
		  const double3& r2,
		  const double3& r3,
		  const double3& r4, // in-plane (ri-rj ^ ri-rk) position
		  const double3& F1,
		  const double3& F2,
		  const double3& F3,
		  double& f_decomp_err,
		  double3* dF) {
  auto dr21 = r1 - r2;
  auto dr41 = r1 - r4;
  auto dr32 = r2 - r3;
  auto dr42 = r2 - r4;
  auto dr43 = r3 - r4;

  constexpr int nRows = 12, nCols = 5, nRhs = 1;

  std::array<double, nRows * nCols> D; D.fill(0.0);
  std::array<double, nRows> b;

  D[nRows * 0 + 0] =  dr21.x; D[nRows * 1 + 0] = dr41.x;
  D[nRows * 0 + 1] =  dr21.y; D[nRows * 1 + 1] = dr41.y;
  D[nRows * 0 + 2] =  dr21.z; D[nRows * 1 + 2] = dr41.z;
  b[0] = F1.x; b[1] = F1.y; b[2] = F1.z;
    
  D[nRows * 0 + 3] = -dr21.x; D[nRows * 2 + 3] = dr32.x; D[nRows * 3 + 3] = dr42.x;
  D[nRows * 0 + 4] = -dr21.y; D[nRows * 2 + 4] = dr32.y; D[nRows * 3 + 4] = dr42.y;
  D[nRows * 0 + 5] = -dr21.z; D[nRows * 2 + 5] = dr32.z; D[nRows * 3 + 5] = dr42.z;
  b[3] = F2.x; b[4] = F2.y; b[5] = F2.z;
    
  D[nRows * 2 + 6] = -dr32.x; D[nRows * 4 + 6] = dr43.x;
  D[nRows * 2 + 7] = -dr32.y; D[nRows * 4 + 7] = dr43.y;
  D[nRows * 2 + 8] = -dr32.z; D[nRows * 4 + 8] = dr43.z;
  b[6] = F3.x; b[7] = F3.y; b[8] = F3.z;
    
  D[nRows * 1 + 9]  = -dr41.x; D[nRows * 3 + 9]  = -dr42.x; D[nRows * 4 + 9]  = -dr43.x;
  D[nRows * 1 + 10] = -dr41.y; D[nRows * 3 + 10] = -dr42.y; D[nRows * 4 + 10] = -dr43.y;
  D[nRows * 1 + 11] = -dr41.z; D[nRows * 3 + 11] = -dr42.z; D[nRows * 4 + 11] = -dr43.z;
  b[9] = 0.0;  b[10] = 0.0; b[11] = 0.0;

  call_dgelsd<nRows, nCols, nRhs>(&D[0], &b[0]);

  const double3 dF21(b[0] * dr21.x, b[0] * dr21.y, b[0] * dr21.z);
  const double3 dF41(b[1] * dr41.x, b[1] * dr41.y, b[1] * dr41.z);
  const double3 dF32(b[2] * dr32.x, b[2] * dr32.y, b[2] * dr32.z);
  const double3 dF42(b[3] * dr42.x, b[3] * dr42.y, b[3] * dr42.z);
  const double3 dF43(b[4] * dr43.x, b[4] * dr43.y, b[4] * dr43.z);

  // check err
  const double3 F1_d = dF21 + dF41;
  const double3 F2_d = dF32 + dF42 - dF21;
  const double3 F3_d = dF43 - dF32;
  const double3 F4_d = -dF41 - dF42 - dF43;
    
  f_decomp_err += (F1_d - F1) * (F1_d - F1);
  f_decomp_err += (F2_d - F2) * (F2_d - F2);
  f_decomp_err += (F3_d - F3) * (F3_d - F3);
  f_decomp_err += F4_d * F4_d;

  dF[0] = dF21;
  dF[1] = dF41;
  dF[2] = dF32;
  dF[3] = dF42;
  dF[4] = dF43;
}

double calc_stress_sum(const double3* dF,
		       const double3* dr,
		       const int num_bond) {
  double sigma_sum = 0.0;
  for (int i = 0; i < num_bond; i++) {
    sigma_sum += dF[i].norm2() * dr[i].norm2();
  }
  return sigma_sum;
}

double calc_stress_sum_cent(const std::array<double, 3> cf,
			    const double3* dr) {
  const double dr_norm[] = {
    dr[0] * dr[0],
    dr[1] * dr[1],
    dr[2] * dr[2]
  };

  const double f_dr_prod[] = {
    cf[0] * dr_norm[0],
    cf[1] * dr_norm[1],
    cf[2] * dr_norm[2]
  };
  
  const auto q = cf[0] * cf[1] + cf[1] * cf[2] + cf[2] * cf[0];

  std::cout << "q value is " << q << std::endl;

  const double ret = ((-cf[1] * cf[0] + cf[2] * cf[0] - cf[1] * cf[2]) * f_dr_prod[0]
		      + (cf[0] * cf[2] - cf[1] * cf[2] - cf[1] * cf[0]) * f_dr_prod[2]
		      + (cf[0] * cf[1] + cf[2] * cf[1] + 3.0 * cf[0] * cf[2]) * f_dr_prod[1]) / std::abs(q);
  return ret;
}

// assume cf[1] < 0.0
double3 get_stress_sum_min_pos(const double3* dr, const double3& base) {
  const auto dr01_norm = dr[0].norm2();
  const auto dr12_norm = dr[1].norm2();

#if 0
  const auto dr02_norm = dr[2].norm2();
  const auto dr01_p_dr12_2 = (dr01_norm + dr12_norm) * (dr01_norm + dr12_norm);
  
  const auto x = -(1.0 + dr12_norm / dr01_norm - std::sqrt((dr01_p_dr12_2 - dr02_norm * dr02_norm) / (dr01_norm * dr01_norm)));
  const auto y = -(1.0 + dr01_norm / dr12_norm - std::sqrt((dr01_p_dr12_2 - dr02_norm * dr02_norm) / (dr12_norm * dr12_norm)));

  // check
  // const auto r2 = dr01_norm * dr01_norm / (dr02_norm * dr02_norm);
  // const auto r3 = dr12_norm * dr12_norm / (dr02_norm * dr02_norm);

  // std::cout << "CHECK\n";
  // std::cout << (1 - y) * r2 * x * x + 2.0 * r2 * x * y + (r3 * y - 1) * y << std::endl;
  // std::cout << (1 - x) * r3 * y * y + 2.0 * r3 * x * y + (r2 * x - 1) * x << std::endl;
  // std::cout << "end check\n";

  const auto a = 1.0 / (1.0 + y / x + y);
  const auto b = 1.0 / (1.0 + x / y + x);

  // std::cout << "coef a b\n";
  // std::cout << a << " " << b << std::endl;

  return base - (a * dr[0] - b * dr[1]);
#else
  const auto dr10_hat = -dr[0] / dr01_norm;
  const auto dr12_hat = dr[1] / dr12_norm;
  auto dr13_hat = (dr10_hat + dr12_hat);
  dr13_hat /= dr13_hat.norm2();
  const double dr13_norm = std::sqrt(dr01_norm * dr12_norm);
  const auto dr13 = dr13_hat * dr13_norm;
  return base + dr13;
#endif
}

void check_congruence_of_two_tri(const double3* pos,
				 const double dist0,
				 const double dist1,
				 const double3& min_pos) {
  const auto a = std::sqrt((min_pos - pos[0]) * (min_pos - pos[0]));
  const auto b = std::sqrt((min_pos - pos[2]) * (min_pos - pos[2]));
  const auto c = std::sqrt(dist0);
  const auto d = std::sqrt((min_pos - pos[1]) * (min_pos - pos[1]));
  const auto e = std::sqrt(dist1);
  
  std::cout << std::setprecision(15);
  std::cout << "|r14/r43|, |r12/r24|, |r24/r23|\n";
  std::cout << a / b << ", " << c / d << ", " << d / e << std::endl;
}
