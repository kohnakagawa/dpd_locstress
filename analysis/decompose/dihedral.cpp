#include "ls_tools.hpp"
#include <random>

namespace {
  constexpr int multip		= 1;
  constexpr double coef_k       = 50.0;  
}

// From lammps/src/dihedral_charmm.cpp
void calc_dihedral(const std::array<double3, 4>& pos,
		   std::array<double3, 4>& force,
		   std::array<double3, 9>& dF,
		   const double3& base_pos,
		   double& stress) {
  // 1st bond
  const auto vb1x = pos[0][0] - pos[1][0], vb1y = pos[0][1] - pos[1][1], vb1z = pos[0][2] - pos[1][2];
  
  // 2nd bond
  const auto vb2x = pos[2][0] - pos[1][0], vb2y = pos[2][1] - pos[1][1], vb2z = pos[2][2] - pos[1][2];
  const auto vb2xm = -vb2x, vb2ym = -vb2y, vb2zm = -vb2z;

  // 3rd bond
  const auto vb3x = pos[3][0] - pos[2][0], vb3y = pos[3][1] - pos[2][1], vb3z = pos[3][2] - pos[2][2];
  
  const auto ax = vb1y * vb2zm - vb1z * vb2ym, ay = vb1z * vb2xm - vb1x * vb2zm, az = vb1x * vb2ym - vb1y * vb2xm;
  const auto bx = vb3y * vb2zm - vb3z * vb2ym, by = vb3z * vb2xm - vb3x * vb2zm, bz = vb3x * vb2ym - vb3y * vb2xm;

  const auto rasq = ax * ax + ay * ay + az * az;
  const auto rbsq = bx * bx + by * by + bz * bz;
  const auto rgsq = vb2xm * vb2xm + vb2ym * vb2ym + vb2zm * vb2zm;
  const auto rg = std::sqrt(rgsq);

  auto rginv = 0.0, ra2inv = 0.0, rb2inv = 0.0;
  if (rg > 0) rginv = 1.0 / rg;
  if (rasq > 0) ra2inv = 1.0 / rasq;
  if (rbsq > 0) rb2inv = 1.0 / rbsq;
  const auto rabinv = std::sqrt(ra2inv * rb2inv);

  auto c = (ax * bx + ay * by + az * bz) * rabinv;
  const double s = rg * rabinv * (ax * vb3x + ay * vb3y + az * vb3z);
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;
  
  const int m = multip;
  double p = 1.0, ddf1 = 0.0, df1 = 0.0;
  for (int i = 0; i < m; i++) {
    ddf1 = p * c - df1 * s;
    df1 = p * s + df1 * c;
    p = ddf1;
  }
  
  // p = p * cos_shift[type] + df1 * sin_shift[type];
  // df1 = df1 * cos_shift[type] - ddf1 * sin_shift[type];
  df1 *= -m;
  p += 1.0;
  
  if (m == 0) {
    // p = 1.0 + cos_shift[type];
    df1 = 0.0;
  }
  
  const auto fg	 = vb1x * vb2xm + vb1y * vb2ym + vb1z * vb2zm;
  const auto hg	 = vb3x * vb2xm + vb3y * vb2ym + vb3z * vb2zm;
  const auto fga = fg * ra2inv * rginv;
  const auto hgb = hg * rb2inv * rginv;
  const auto gaa = -ra2inv * rg;
  const auto gbb = rb2inv * rg;

  const auto dtfx = gaa * ax, dtfy = gaa * ay, dtfz = gaa * az;
  const auto dtgx = fga * ax - hgb * bx, dtgy = fga * ay - hgb * by, dtgz = fga * az - hgb *bz;
  const auto dthx = gbb * bx, dthy = gbb * by, dthz = gbb * bz;

  const auto df = -coef_k * df1;

  const auto sx2 = df * dtgx;
  const auto sy2 = df * dtgy;
  const auto sz2 = df * dtgz;

  const double3 f1(df * dtfx, df * dtfy, df * dtfz);
  const double3 f2(sx2 - f1[0], sy2 - f1[1], sz2 - f1[2]);
  const double3 f4(df * dthx, df * dthy, df * dthz);
  const double3 f3(-sx2 - f4[0], -sy2 - f4[1], -sz2 - f4[2]);

  std::cout << "linear momentum.\n";
  std::cout << f1 + f2 + f3 + f4 << std::endl;
  std::cout << "angular momentum.\n";
  std::cout << (f1 ^ pos[0]) + (f2 ^ pos[1]) + (f3 ^ pos[2]) + (f4 ^ pos[3]) << std::endl;

  double f_decomp_err = 0.0;
  decompose_4n(pos[0], pos[1], pos[3], pos[2], base_pos, f1, f2, f4, f3, f_decomp_err, stress, &dF[0]);

  std::cout << "decompose error\n";
  std::cout << f_decomp_err << std::endl;
  
  // apply force to each of 4 atoms
  force[0] += f1;
  force[1] += f2;
  force[2] += f3;
  force[3] += f4;
}

void gen_init_config(std::array<double3, 4>& pos,
		     const std::array<double, 3>& b_len,
		     const double angle,
		     const double phi) {
  pos[0] = 0.0;
  
  pos[1][0] = b_len[0] * std::sin(angle);
  pos[1][1] = 0.0;
  pos[1][2] = b_len[0] * std::cos(angle);

  pos[2][0] = pos[1][0];
  pos[2][1] = pos[1][1];
  pos[2][2] = pos[1][2] + b_len[1];

  pos[3][0] = pos[2][0] + b_len[2] * std::sin(angle) * std::cos(phi);
  pos[3][1] = pos[2][1] + b_len[2] * std::sin(angle) * std::sin(phi);
  pos[3][2] = pos[2][2] + b_len[2] * std::cos(angle);

  const double3 shift(5.0);

  for (auto i = 0u; i < pos.size(); i++)
    pos[i] += shift;

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      assert(pos[i][j] >= 0.0);

  std::cout << "particle configuration.\n";
  for (int i = 0; i < 4; i++)
    std::cout << pos[i] << std::endl;
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "argv[1] == target directory name.\n";
    std::cerr << "argv[2] == dihedral angle.\n";
    std::exit(1);
  }
  const double phi = std::atof(argv[2]);

  std::array<double3, 4> pos, force;
  std::array<double3, 9> dF;
  std::array<double, 3> b_len{1.0, 1.0, 1.0};
  const double angle = M_PI * 30.0 / 180.0;
  
  double stress = 0.0;
  gen_init_config(pos, b_len, angle, phi);
  // const double3 base_pos = (pos[0] + pos[2]) * 0.5;
  const double3 base_pos = std::accumulate(pos.cbegin(), pos.cend(), double3(0.0)) / pos.size();
  force.fill(double3(0.0));
  dF.fill(double3(0.0));
  calc_dihedral(pos, force, dF, base_pos, stress);
}
