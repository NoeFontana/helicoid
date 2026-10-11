// Hand-computed tests of every conversion in `convert.hpp` and of the hex floats, on values whose
// components all differ, so a permuted, transposed or reversed layout fails here rather than scoring
// (docs/PHASE1.md §7). The image build runs this and stops on a failure.
#include <cmath>
#include <cstdio>
#include <cstring>
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

// Far above the rounding of the inputs and of Sophus at these values, far below any slip.
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

void quaternions() {
  const Eigen::Quaterniond q = oracle::helicoid_to_sophus_quat({0.1, 0.2, 0.3, 0.4});
  check(q.w() == 0.1 && q.x() == 0.2 && q.y() == 0.3 && q.z() == 0.4, "quat: (w, x, y, z)");
  check(q.coeffs()(3) == 0.1 && q.coeffs()(0) == 0.2, "quat: Eigen stores w last");
  exact(oracle::sophus_to_helicoid_quat(q), {0.1, 0.2, 0.3, 0.4}, "quat: back");
}

void tangents() {
  using oracle::helicoid_to_sophus_tangent;
  using oracle::sophus_to_helicoid_tangent;
  const auto se2 = helicoid_to_sophus_tangent({1, 2, 3}, oracle::kSe2);
  check(se2 == Eigen::Vector3d(2, 3, 1), "se2 tangent: [θ; ρ] → [υ; θ]");
  exact(sophus_to_helicoid_tangent(se2, oracle::kSe2), {1, 2, 3}, "se2 tangent: back");
  const auto se3 = helicoid_to_sophus_tangent({1, 2, 3, 4, 5, 6}, oracle::kSe3);
  check(se3 == (Sophus::SE3d::Tangent() << 4, 5, 6, 1, 2, 3).finished(), "se3 tangent: [φ; ρ] → [υ; ω]");
  exact(sophus_to_helicoid_tangent(se3, oracle::kSe3), {1, 2, 3, 4, 5, 6}, "se3 tangent: back");
  const auto sim3 = helicoid_to_sophus_tangent({1, 2, 3, 4, 5, 6, 7}, oracle::kSim3);
  check(sim3 == (Sophus::Sim3d::Tangent() << 4, 5, 6, 1, 2, 3, 7).finished(),
        "sim3 tangent: [φ; ρ; σ] → [υ; ω; σ]");
  exact(sophus_to_helicoid_tangent(sim3, oracle::kSim3), {1, 2, 3, 4, 5, 6, 7}, "sim3 tangent: back");
}

void matrices() {
  // m(r, c) = 10 r + c in Sophus's order [υ₀; υ₁; θ]; helicoid's (i, j) is Sophus's (at[i], at[j]).
  Eigen::Matrix3d m;
  m << 0, 1, 2, 10, 11, 12, 20, 21, 22;
  exact(oracle::sophus_to_helicoid_matrix(m, oracle::kSe2), {22, 2, 12, 20, 0, 10, 21, 1, 11},
        "se2 matrix: permuted, column-major");
  exact(oracle::sophus_to_helicoid_matrix3(m), {0, 10, 20, 1, 11, 21, 2, 12, 22}, "matrix3: column-major");
}

void se2_adjoint() {
  // θ = π/2, t = (1, 2): helicoid's Ad = [1 0 0; t₁ c −s; −t₀ s c] = [1 0 0; 2 0 −1; −1 1 0].
  const Sophus::SE2d g = oracle::helicoid_to_sophus_se2({0, 1}, {1, 2});
  near(oracle::sophus_to_helicoid_matrix(g.Adj(), oracle::kSe2), {1, 2, -1, 0, 0, 1, 0, -1, 0}, "se2 Ad");
  near(oracle::sophus_to_helicoid_so2(oracle::helicoid_to_sophus_so2({0.6, 0.8})), {0.6, 0.8}, "so2: (re, im)");
}

void se3_left_jacobian() {
  // φ = 0, ρ = (1, 2, 3): helicoid's J_l = [I 0; ½ρ^ I], ½ρ^ = [0 −1.5 1; 1.5 0 −0.5; −1 0.5 0].
  const auto tau = oracle::helicoid_to_sophus_tangent({0, 0, 0, 1, 2, 3}, oracle::kSe3);
  exact(oracle::sophus_to_helicoid_matrix(Sophus::SE3d::leftJacobian(tau), oracle::kSe3),
        {1, 0, 0, 0, 1.5, -1, 0, 1, 0, -1.5, 0, 0.5, 0, 0, 1, 1, -0.5, 0,
         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        "se3 J_l at φ = 0");
}

void se3_element() {
  const double h = std::sqrt(0.5);
  const Sophus::SE3d g = oracle::helicoid_to_sophus_se3({h, 0, 0, h}, {1, 2, 3});
  near(oracle::sophus_to_helicoid_vec3(g * Eigen::Vector3d(1, 0, 0)), {1, 3, 3}, "se3: R_z(π/2) e_x + t");
}

void sim3_element() {
  // q = (1, 2, 3, 4)/√30, σ = ln 4: Sophus stores 2q, and reads σ back from its squared norm.
  const double r = std::sqrt(30.0);
  const oracle::Vec q{1 / r, 2 / r, 3 / r, 4 / r};
  const Sophus::Sim3d x = oracle::helicoid_to_sophus_sim3(q, {5, 6, 7}, std::log(4.0));
  const Eigen::Quaterniond& stored = x.rxso3().quaternion();
  near({stored.w(), stored.x(), stored.y(), stored.z()}, {2 / r, 4 / r, 6 / r, 8 / r}, "sim3: √s · q");
  near(oracle::sophus_to_helicoid_vec3(x.translation()), {5, 6, 7}, "sim3: x");
  const oracle::HelicoidSim3 back = oracle::sophus_to_helicoid_sim3(x);
  near(back.q, q, "sim3: q back");
  near({back.sigma}, {std::log(4.0)}, "sim3: σ back");
  bool refused = false;
  try {
    oracle::helicoid_to_sophus_sim3(q, {5, 6, 7}, 30.0);
  } catch (const oracle::Refused&) {
    refused = true;
  }
  check(refused, "sim3: e^30 is past RxSO3's 10^10, refused");
}

void sim3_adjoint() {
  // R = R_z(π/2), s = 2, t = (1, 2, 3): helicoid's Ad = [R 0 0; t^R sR −t; 0 0 1], with
  // t^R = [−3 0 2; 0 −3 −1; 1 2 0] and sR = [0 −2 0; 2 0 0; 0 0 2].
  const double h = std::sqrt(0.5);
  const Sophus::Sim3d x = oracle::helicoid_to_sophus_sim3({h, 0, 0, h}, {1, 2, 3}, std::log(2.0));
  near(oracle::sophus_to_helicoid_matrix(x.Adj(), oracle::kSim3),
       {0, 1, 0, -3, 0, 1, 0, -1, 0, 0, 0, -3, 2, 0, 0, 0, 1, 2, -1, 0, 0, 0, 0, 0, 0, 2, 0, 0,
        0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, -1, -2, -3, 1},
       "sim3 Ad");
}

void sim3_exp() {
  // [φ; ρ; σ] = [(0, 0, π/2); 0; ln 2] is (q, x, σ) = ((√½, 0, 0, √½), 0, ln 2); with φ = σ = 0,
  // x = ρ.
  const double h = std::sqrt(0.5);
  const auto turn = oracle::sophus_to_helicoid_sim3(Sophus::Sim3d::exp(
      oracle::helicoid_to_sophus_tangent({0, 0, 2 * std::atan(1.0), 0, 0, 0, std::log(2.0)}, oracle::kSim3)));
  near(turn.q, {h, 0, 0, h}, "sim3 exp: q");
  near(turn.x, {0, 0, 0}, "sim3 exp: x");
  near({turn.sigma}, {std::log(2.0)}, "sim3 exp: σ");
  const auto shift = oracle::sophus_to_helicoid_sim3(
      Sophus::Sim3d::exp(oracle::helicoid_to_sophus_tangent({0, 0, 0, 1, 2, 3, 0}, oracle::kSim3)));
  near(shift.x, {1, 2, 3}, "sim3 exp: x = ρ");
}

void hex_floats() {
  check(oracle::format_hex(1.0) == "0x1p+0", "hex: 1");
  check(oracle::format_hex(-0.0) == "-0x0p+0", "hex: −0");
  check(oracle::format_hex(4.9406564584124654e-324) == "0x0.0000000000001p-1022", "hex: subnormal");
  check(oracle::format_hex(std::nan("")) == "nan" && oracle::format_hex(-INFINITY) == "-inf", "hex: words");
  check(same_bits(oracle::parse_hex("0x1.921fb54442d18p+0"), 2 * std::atan(1.0)), "hex: π/2");
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
  se2_adjoint();
  se3_left_jacobian();
  se3_element();
  sim3_element();
  sim3_adjoint();
  sim3_exp();
  hex_floats();
  if (failures != 0) return 1;
  std::puts("convert_test: all passed");
  return 0;
}
