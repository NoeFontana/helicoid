// Hand-computed tests of every conversion in `convert.hpp` and of the hex floats, on values whose
// components all differ, so a permuted, transposed or reversed layout fails here rather than scoring
// (docs/PHASE1.md §7). The image build runs this and stops on a failure.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <initializer_list>
#include <string>

#include "convert.hpp"
#include "hexfloat.hpp"

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", what.c_str());
    ++failures;
  }
}

// Far above the rounding of the inputs and of manif at these values, far below any slip.
constexpr double kTol = 1e-14;

void near(const oracle::Vec& got, const oracle::Vec& want, const std::string& what) {
  bool ok = got.size() == want.size();
  for (std::size_t i = 0; ok && i < got.size(); ++i) ok = std::fabs(got[i] - want[i]) <= kTol;
  check(ok, what);
}

void exact(const oracle::Vec& got, const oracle::Vec& want, const std::string& what) {
  check(got == want, what);
}

bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

const double kH = std::sqrt(0.5);
const double kHalfPi = 2 * std::atan(1.0);

// A rotation-first matrix of 3 × 3 blocks, each written row by row, as the corpus's array:
// column-major. `blocks[r][c]` is block row `r`, block column `c`.
oracle::Vec blocks3(int n, std::initializer_list<std::initializer_list<oracle::Vec>> blocks) {
  oracle::Vec out(9 * n * n, 0.0);
  int r = 0;
  for (const auto& row : blocks) {
    int c = 0;
    for (const auto& b : row) {
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) out[(3 * c + j) * 3 * n + 3 * r + i] = b.at(3 * i + j);
      }
      ++c;
    }
    ++r;
  }
  return out;
}

const oracle::Vec kI{1, 0, 0, 0, 1, 0, 0, 0, 1};
const oracle::Vec k0(9, 0.0);
// R_z(π/2).
const oracle::Vec kRz{0, -1, 0, 1, 0, 0, 0, 0, 1};

void quaternions() {
  const Eigen::Quaterniond q = oracle::helicoid_to_manif_quat({0.1, 0.2, 0.3, 0.4});
  check(q.w() == 0.1 && q.x() == 0.2 && q.y() == 0.3 && q.z() == 0.4, "quat: (w, x, y, z)");
  check(q.coeffs()(3) == 0.1 && q.coeffs()(0) == 0.2, "quat: Eigen stores w last");
  exact(oracle::manif_to_helicoid_quat(q), {0.1, 0.2, 0.3, 0.4}, "quat: back");
  const manif::SO3d r = oracle::helicoid_to_manif_so3({kH, 0, 0, kH});
  check(r.coeffs()(3) == kH && r.coeffs()(2) == kH && r.coeffs()(0) == 0, "so3: manif stores [x, y, z, w]");
}

void tangents() {
  using oracle::helicoid_to_manif_tangent;
  using oracle::manif_to_helicoid_tangent;
  const auto se2 = helicoid_to_manif_tangent({1, 2, 3}, oracle::kSe2);
  check(se2 == Eigen::Vector3d(2, 3, 1), "se2 tangent: [θ; ρ] → [ρ; θ]");
  exact(manif_to_helicoid_tangent(se2, oracle::kSe2), {1, 2, 3}, "se2 tangent: back");
  const auto se3 = helicoid_to_manif_tangent({1, 2, 3, 4, 5, 6}, oracle::kSe3);
  check(se3 == (manif::SE3d::Tangent::DataType() << 4, 5, 6, 1, 2, 3).finished(), "se3 tangent: [φ; ρ] → [ρ; θ]");
  exact(manif_to_helicoid_tangent(se3, oracle::kSe3), {1, 2, 3, 4, 5, 6}, "se3 tangent: back");
  const auto se23 = helicoid_to_manif_tangent({1, 2, 3, 4, 5, 6, 7, 8, 9}, oracle::kSe23);
  check(se23 == (manif::SE_2_3d::Tangent::DataType() << 7, 8, 9, 1, 2, 3, 4, 5, 6).finished(),
        "se23 tangent: [φ; ρ_v; ρ_p] → [ρ_p; θ; ν]");
  exact(manif_to_helicoid_tangent(se23, oracle::kSe23), {1, 2, 3, 4, 5, 6, 7, 8, 9}, "se23 tangent: back");
}

void matrices() {
  // m(r, c) = 10 r + c in manif's order [ρ₀; ρ₁; θ]; helicoid's (i, j) is manif's (at[i], at[j]).
  Eigen::Matrix3d m;
  m << 0, 1, 2, 10, 11, 12, 20, 21, 22;
  exact(oracle::manif_to_helicoid_matrix(m, oracle::kSe2), {22, 2, 12, 20, 0, 10, 21, 1, 11},
        "se2 matrix: permuted, column-major");
  exact(oracle::manif_to_helicoid_matrix3(m), {0, 10, 20, 1, 11, 21, 2, 12, 22}, "matrix3: column-major");
  exact(blocks3(1, {{{0, 1, 2, 10, 11, 12, 20, 21, 22}}}), {0, 10, 20, 1, 11, 21, 2, 12, 22}, "blocks3: column-major");
}

void so2_and_se2() {
  near(oracle::manif_to_helicoid_so2(oracle::helicoid_to_manif_so2({0.6, 0.8})), {0.6, 0.8}, "so2: (re, im)");
  near({oracle::helicoid_to_manif_so2({0.6, 0.8}).angle()}, {std::atan2(0.8, 0.6)}, "so2: angle of (re, im)");
  bool refused = false;
  try {
    oracle::helicoid_to_manif_so2({0.6, 0.81});
  } catch (const manif::invalid_argument&) {
    refused = true;
  }
  check(refused, "so2: a non-unit complex is refused, not normalized");
  // θ = π/2, t = (1, 2): helicoid's Ad = [1 0 0; t₁ c −s; −t₀ s c] = [1 0 0; 2 0 −1; −1 1 0].
  const manif::SE2d g = oracle::helicoid_to_manif_se2({0, 1}, {1, 2});
  near(oracle::manif_to_helicoid_matrix(g.adj(), oracle::kSe2), {1, 2, -1, 0, 0, 1, 0, -1, 0}, "se2 Ad");
  const oracle::HelicoidSe2 back = oracle::manif_to_helicoid_se2(g);
  near(back.z, {0, 1}, "se2: z back");
  near(back.t, {1, 2}, "se2: t back");
  // θ = 0, ρ = (1, 2): helicoid's J_l = I + ½ ad_τ = [1 0 0; ρ₁/2 1 0; −ρ₀/2 0 1], J_r = I − ½ ad_τ.
  const auto tau = oracle::helicoid_to_manif_tangent({0, 1, 2}, oracle::kSe2);
  const manif::SE2Tangentd t(tau);
  near(oracle::manif_to_helicoid_matrix(t.ljac(), oracle::kSe2), {1, 1, -0.5, 0, 1, 0, 0, 0, 1}, "se2 J_l at θ = 0");
  near(oracle::manif_to_helicoid_matrix(t.rjac(), oracle::kSe2), {1, -1, 0.5, 0, 1, 0, 0, 0, 1}, "se2 J_r at θ = 0");
}

void so3() {
  // φ = (0, 0, π/2): J_l = I + (1 − cos θ)/θ² W + (θ − sin θ)/θ³ W² = [2/π −2/π 0; 2/π 2/π 0; 0 0 1],
  // and J_r is its transpose.
  const double a = 1 / kHalfPi;
  const manif::SO3Tangentd phi(Eigen::Vector3d(0, 0, kHalfPi));
  near(oracle::manif_to_helicoid_matrix3(phi.ljac()), {a, a, 0, -a, a, 0, 0, 0, 1}, "so3 J_l");
  near(oracle::manif_to_helicoid_matrix3(phi.rjac()), {a, -a, 0, a, a, 0, 0, 0, 1}, "so3 J_r");
  near(oracle::manif_to_helicoid_quat(phi.exp().quat()), {kH, 0, 0, kH}, "so3 exp: q = (cos π/4, 0, 0, sin π/4)");
  const manif::SO3d r = oracle::helicoid_to_manif_so3({kH, 0, 0, kH});
  near(oracle::manif_to_helicoid_vec3(r.act(Eigen::Vector3d(1, 0, 0))), {0, 1, 0}, "so3 act: R_z(π/2) e_x");
  near(oracle::manif_to_helicoid_vec3(r.log().coeffs()), {0, 0, kHalfPi}, "so3 log");
  const manif::SO3d mid = manif::interpolate_slerp(manif::SO3d::Identity(), r, 0.5);
  near(oracle::manif_to_helicoid_quat(mid.quat()), {std::cos(kHalfPi / 4), 0, 0, std::sin(kHalfPi / 4)},
       "so3 geodesic: half of R_z(π/2)");
  bool refused = false;
  try {
    oracle::helicoid_to_manif_so3({1, 1e-6, 0, 0});
  } catch (const manif::invalid_argument&) {
    refused = true;
  }
  check(refused, "so3: ‖q‖ − 1 = 5·10⁻¹³, past 100 ε, is refused, not normalized");
}

void se3() {
  const manif::SE3d g = oracle::helicoid_to_manif_se3({kH, 0, 0, kH}, {1, 2, 3});
  near(oracle::manif_to_helicoid_vec3(g.act(Eigen::Vector3d(1, 0, 0))), {1, 3, 3}, "se3: R_z(π/2) e_x + t");
  const oracle::HelicoidSen3 back = oracle::manif_to_helicoid_se3(g);
  near(back.q, {kH, 0, 0, kH}, "se3: q back");
  near(back.x, {1, 2, 3}, "se3: x back");
  // φ = 0, ρ = (1, 2, 3): helicoid's J_l = [I 0; ½ρ^ I], ½ρ^ = [0 −1.5 1; 1.5 0 −0.5; −1 0.5 0].
  const manif::SE3Tangentd tau(oracle::helicoid_to_manif_tangent({0, 0, 0, 1, 2, 3}, oracle::kSe3));
  const oracle::Vec half{0, -1.5, 1, 1.5, 0, -0.5, -1, 0.5, 0};
  exact(oracle::manif_to_helicoid_matrix(tau.ljac(), oracle::kSe3), blocks3(2, {{kI, k0}, {half, kI}}),
        "se3 J_l at φ = 0");
  // helicoid's Ad = [R 0; [t]×R R], [t]×R = [−3 0 2; 0 −3 −1; 1 2 0].
  near(oracle::manif_to_helicoid_matrix(g.adj(), oracle::kSe3),
       blocks3(2, {{kRz, k0}, {{-3, 0, 2, 0, -3, -1, 1, 2, 0}, kRz}}), "se3 Ad");
  // The screw chart: X Exp([0; (1, 0, 0)]) = (R, t + R e_x).
  const manif::SE3Tangentd shift(oracle::helicoid_to_manif_tangent({0, 0, 0, 1, 0, 0}, oracle::kSe3));
  near(oracle::manif_to_helicoid_se3(g.rplus(shift)).x, {1, 3, 3}, "se3 screw retract: t + R ρ");
  near(oracle::manif_to_helicoid_tangent(g.rplus(shift).rminus(g).coeffs(), oracle::kSe3), {0, 0, 0, 1, 0, 0},
       "se3 screw local");
}

void world_chart() {
  // `(q Exp φ, x + ρ)`: the translation is not rotated, unlike the screw and decoupled charts.
  const oracle::SO3xR3 x = oracle::helicoid_to_manif_so3xr3({kH, 0, 0, kH}, {1, 2, 3});
  const oracle::SO3xR3::Tangent tau((oracle::SO3xR3::Tangent::DataType() << 0, 0, kHalfPi, 1, 0, 0).finished());
  const oracle::HelicoidSen3 y = oracle::manif_to_helicoid_so3xr3(x.rplus(tau));
  near(y.q, {0, 0, 0, 1}, "world retract: q Exp φ = R_z(π)");
  near(y.x, {2, 2, 3}, "world retract: x + ρ");
  const oracle::SO3xR3 z = oracle::helicoid_to_manif_so3xr3({kH, 0, 0, kH}, {4, 6, 8});
  const Eigen::VectorXd back = z.rminus(x).coeffs();
  near(oracle::Vec(back.data(), back.data() + 6), {0, 0, 0, 3, 4, 5}, "world local: [Log(q₀* q₁); x₁ − x₀]");
}

void se23() {
  // helicoid's (q, [v; p]) with v = (1, 2, 3), p = (4, 5, 6): manif's t is p, its v is v.
  const manif::SE_2_3d g = oracle::helicoid_to_manif_se23({kH, 0, 0, kH}, {1, 2, 3, 4, 5, 6});
  check(g.translation() == Eigen::Vector3d(4, 5, 6), "se23: manif's translation is helicoid's p = x₂");
  check(g.linearVelocity() == Eigen::Vector3d(1, 2, 3), "se23: manif's velocity is helicoid's v = x₁");
  const oracle::HelicoidSen3 back = oracle::manif_to_helicoid_se23(g);
  near(back.q, {kH, 0, 0, kH}, "se23: q back");
  exact(back.x, {1, 2, 3, 4, 5, 6}, "se23: x back");
  // helicoid's Ad = [R 0 0; [v]×R R 0; [p]×R 0 R] (docs/NUMERICS.md §5.2), [v]×R = [−3 0 2; 0 −3 −1;
  // 1 2 0], [p]×R = [−6 0 5; 0 −6 −4; 4 5 0].
  near(oracle::manif_to_helicoid_matrix(g.adj(), oracle::kSe23),
       blocks3(3, {{kRz, k0, k0}, {{-3, 0, 2, 0, -3, -1, 1, 2, 0}, kRz, k0}, {{-6, 0, 5, 0, -6, -4, 4, 5, 0}, k0, kRz}}),
       "se23 Ad");
  // φ = 0, ρ_v = (1, 2, 3), ρ_p = (4, 5, 6): J_l = [I 0 0; ½ρ_v^ I 0; ½ρ_p^ 0 I], and Exp is
  // (I, ρ_v, ρ_p).
  const manif::SE_2_3Tangentd tau(oracle::helicoid_to_manif_tangent({0, 0, 0, 1, 2, 3, 4, 5, 6}, oracle::kSe23));
  exact(oracle::manif_to_helicoid_matrix(tau.ljac(), oracle::kSe23),
        blocks3(3, {{kI, k0, k0}, {{0, -1.5, 1, 1.5, 0, -0.5, -1, 0.5, 0}, kI, k0}, {{0, -3, 2.5, 3, 0, -2, -2.5, 2, 0}, k0, kI}}),
        "se23 J_l at φ = 0");
  exact(oracle::manif_to_helicoid_se23(tau.exp()).x, {1, 2, 3, 4, 5, 6}, "se23 exp at φ = 0: x = [ρ_v; ρ_p]");
}

void hex_floats() {
  check(oracle::format_hex(1.0) == "0x1p+0", "hex: 1");
  check(oracle::format_hex(-0.0) == "-0x0p+0", "hex: −0");
  check(oracle::format_hex(4.9406564584124654e-324) == "0x0.0000000000001p-1022", "hex: subnormal");
  check(oracle::format_hex(std::nan("")) == "nan" && oracle::format_hex(-INFINITY) == "-inf", "hex: words");
  check(same_bits(oracle::parse_hex("0x1.921fb54442d18p+0"), kHalfPi), "hex: π/2");
  check(same_bits(oracle::parse_hex("-0x0.0p+0"), -0.0), "hex: −0 in");
  check(oracle::parse_hex("0x0.0000000000001p-1022") == 4.9406564584124654e-324, "hex: subnormal in");
  for (const char* bad : {"1.5", "0x1p+0 ", "nan", "0x1.8", ""}) {
    bool threw = false;
    try {
      oracle::parse_hex(bad);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    check(threw, std::string("hex: refuses `") + bad + "`");
  }
}

}  // namespace

int main() {
  quaternions();
  tangents();
  matrices();
  so2_and_se2();
  so3();
  se3();
  world_chart();
  se23();
  hex_floats();
  if (failures != 0) return 1;
  std::puts("convert_test: all passed");
  return 0;
}
