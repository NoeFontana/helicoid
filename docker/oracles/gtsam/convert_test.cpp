// Hand-computed tests of every conversion in `convert.hpp`, of the side of each GTSAM derivative
// the runner reads, and of the hex floats, on values whose components all differ, so a permuted,
// transposed, reversed or wrong-sided layout fails here rather than scoring (docs/PHASE1.md §7). The
// image build runs this and stops on a failure.
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

// Far above the rounding of the inputs and of GTSAM at these values, far below any slip.
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

const double kHalf = std::sqrt(0.5);
const double kPi = 4 * std::atan(1.0);

void rotations() {
  // q = (1, 2, 3, 4)/√30: R(0, 1) = 2(xy − wz) = 4/30, R(1, 0) = 2(xy + wz) = 20/30.
  const double r = std::sqrt(30.0);
  const gtsam::Rot3 g = oracle::helicoid_to_gtsam_rot3({1 / r, 2 / r, 3 / r, 4 / r});
  check(std::fabs(g.matrix()(0, 1) - 4.0 / 30) <= kTol && std::fabs(g.matrix()(1, 0) - 20.0 / 30) <= kTol,
        "rot3: (w, x, y, z), Hamilton, active");
  near(oracle::gtsam_to_helicoid_quat(g), {1 / r, 2 / r, 3 / r, 4 / r}, "rot3: quaternion back, w first");
  // q = (0, 2, 0, 0) is not normalized: `toRotationMatrix` gives diag(1, 1 − 8, 1 − 8).
  const gtsam::Rot3 unnormalized = oracle::helicoid_to_gtsam_rot3({0, 2, 0, 0});
  check(unnormalized.matrix()(1, 1) == -7 && unnormalized.matrix()(0, 0) == 1, "rot3: not normalized");
  near(oracle::gtsam_to_helicoid_vec3(oracle::helicoid_to_gtsam_rot3({kHalf, 0, 0, kHalf}).rotate(Eigen::Vector3d(1, 0, 0))),
       {0, 1, 0}, "rot3: R_z(π/2) e_x = e_y");
  // The corpus's `R` is column-major.
  const Eigen::Matrix3d m = oracle::helicoid_to_gtsam_matrix3({0, 1, 2, 3, 4, 5, 6, 7, 8});
  check(m(1, 0) == 1 && m(0, 1) == 3 && m(2, 2) == 8, "matrix3: column-major in");
}

void tangents() {
  const Eigen::Vector3d se2 = oracle::helicoid_to_gtsam_se2_tangent({1, 2, 3});
  check(se2 == Eigen::Vector3d(2, 3, 1), "se2 tangent: [θ; ρ] → [v; ω]");
  exact(oracle::gtsam_to_helicoid_se2_tangent(se2), {1, 2, 3}, "se2 tangent: back");
  const auto sim3 = oracle::helicoid_to_gtsam_tangent<7>({1, 2, 3, 4, 5, 6, 7});
  check(sim3 == (gtsam::Vector7() << 1, 2, 3, 4, 5, 6, 7).finished(), "sim3 tangent: [φ; ρ; σ] = [ω; ρ; λ]");
  exact(oracle::gtsam_to_helicoid_tangent<7>(sim3), {1, 2, 3, 4, 5, 6, 7}, "sim3 tangent: back");
}

void matrices() {
  // m(r, c) = 10 r + c in GTSAM's order [v₀; v₁; ω]; helicoid's (i, j) is GTSAM's (at[i], at[j]).
  Eigen::Matrix3d m;
  m << 0, 1, 2, 10, 11, 12, 20, 21, 22;
  exact(oracle::gtsam_to_helicoid_se2_matrix(m), {22, 2, 12, 20, 0, 10, 21, 1, 11},
        "se2 matrix: permuted, column-major");
  exact(oracle::gtsam_to_helicoid_matrix(m), {0, 10, 20, 1, 11, 21, 2, 12, 22}, "matrix: column-major");
}

void se2() {
  // θ = π/2, t = (1, 2): helicoid's Ad = [1 0 0; t₁ c −s; −t₀ s c] = [1 0 0; 2 0 −1; −1 1 0].
  const gtsam::Pose2 g = oracle::helicoid_to_gtsam_pose2({0, 1}, {1, 2});
  near(oracle::gtsam_to_helicoid_se2_matrix(g.AdjointMap()), {1, 2, -1, 0, 0, 1, 0, -1, 0}, "se2 Ad");
  near(oracle::gtsam_to_helicoid_rot2(oracle::helicoid_to_gtsam_rot2({0.6, 0.8})), {0.6, 0.8}, "so2: (c, s)");
}

void so3_sides() {
  // φ = (0, 0, π/2): J_r = I − aW + bW², a = 4/π², so J_r(0, 1) = −a W(0, 1) = 2/π, and J_l's
  // is −2/π. J_r⁻¹ = I + ½W + cW², so J_r⁻¹(0, 1) = −π/4.
  const Eigen::Vector3d phi(0, 0, kPi / 2);
  check(std::fabs(gtsam::Rot3::ExpmapDerivative(phi)(0, 1) - 2 / kPi) <= kTol, "Rot3::ExpmapDerivative is J_r");
  check(std::fabs(gtsam::Rot3::LogmapDerivative(phi)(0, 1) + kPi / 4) <= kTol, "Rot3::LogmapDerivative is J_r⁻¹");
  // Γ₂ = ½I + bW + dW²: at the same φ, Γ₂(0, 1) = −b π/2 with b = (θ − sin θ)/θ³.
  const gtsam::so3::DexpFunctor local(phi);
  const double b = (kPi / 2 - 1) / std::pow(kPi / 2, 3);
  check(std::fabs(local.Gamma().left()(0, 1) + b * kPi / 2) <= kTol, "DexpFunctor::Gamma().left() is Γ₂");
}

void se3() {
  const gtsam::Pose3 g = oracle::helicoid_to_gtsam_pose3({kHalf, 0, 0, kHalf}, {1, 2, 3});
  near(oracle::gtsam_to_helicoid_vec3(g.transformFrom(Eigen::Vector3d(1, 0, 0))), {1, 3, 3}, "se3: R_z(π/2) e_x + t");
  // φ = 0, ρ = (1, 2, 3): helicoid's J_r = [I 0; −½ρ^ I], −½ρ^ = [0 1.5 −1; −1.5 0 0.5; 1 −0.5 0],
  // and J_r⁻¹ = [I 0; ½ρ^ I].
  const gtsam::Vector6 tau = oracle::helicoid_to_gtsam_tangent<6>({0, 0, 0, 1, 2, 3});
  near(oracle::gtsam_to_helicoid_matrix(gtsam::Pose3::ExpmapDerivative(tau)),
        {1, 0, 0, 0, -1.5, 1, 0, 1, 0, 1.5, 0, -0.5, 0, 0, 1, -1, 0.5, 0,
         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        "Pose3::ExpmapDerivative is J_r, rotation-first");
  near(oracle::gtsam_to_helicoid_matrix(gtsam::Pose3::LogmapDerivative(tau)),
        {1, 0, 0, 0, 1.5, -1, 0, 1, 0, -1.5, 0, 0.5, 0, 0, 1, 1, -0.5, 0,
         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        "Pose3::LogmapDerivative is J_r⁻¹");
}

void extended_pose3() {
  // x = [x₁; x₂] = [(1, 2, 3); (4, 5, 6)] at R = I: the columns in order, and Ad's lower-left
  // blocks [x₁]× and [x₂]×, so Ad(4, 0) = x₁z = 3 and Ad(7, 0) = x₂z = 6.
  const gtsam::ExtendedPose3<2> g = oracle::helicoid_to_gtsam_extended_pose3<2>({1, 0, 0, 0}, {1, 2, 3, 4, 5, 6});
  exact(oracle::gtsam_to_helicoid_matrix(g.xMatrix()), {1, 2, 3, 4, 5, 6}, "se23: columns in order");
  const auto ad = g.AdjointMap();
  check(ad(4, 0) == 3 && ad(7, 0) == 6 && ad(3, 1) == -3 && ad(6, 1) == -6, "se23 Ad: [x_i]× R");
  // J_r at φ = 0 has −½[ρ_i]× in block row i: ρ₂ = (4, 5, 6), so J_r(7, 0) = −½ρ₂z = −3.
  const auto j = gtsam::ExtendedPose3<2>::ExpmapDerivative(
      oracle::helicoid_to_gtsam_tangent<9>({0, 0, 0, 1, 2, 3, 4, 5, 6}));
  near({j(7, 0), j(4, 0), j(7, 6), j(6, 6)}, {-3, -1.5, 0, 1}, "ExtendedPose3::ExpmapDerivative is J_r");
}

void sim3_element() {
  // q = (1, 2, 3, 4)/√30, x = (5, 6, 7), σ = ln 4: GTSAM stores t = x/4 and s = 4.
  const double r = std::sqrt(30.0);
  const oracle::Vec q{1 / r, 2 / r, 3 / r, 4 / r};
  const gtsam::Similarity3 g = oracle::helicoid_to_gtsam_similarity3(q, {5, 6, 7}, std::log(4.0));
  near(oracle::gtsam_to_helicoid_vec3(g.translation()), {1.25, 1.5, 1.75}, "sim3: t = x/s");
  near({g.scale()}, {4}, "sim3: s = e^σ");
  const oracle::HelicoidSim3 back = oracle::gtsam_to_helicoid_sim3(g);
  near(back.q, q, "sim3: q back");
  near(back.x, {5, 6, 7}, "sim3: x = s t back");
  near({back.sigma}, {std::log(4.0)}, "sim3: σ back");
  // R = R_z(π/2), s = 2, x = (1, 2, 3): X e_x = s R e_x + x = (1, 4, 3).
  const gtsam::Similarity3 h = oracle::helicoid_to_gtsam_similarity3({kHalf, 0, 0, kHalf}, {1, 2, 3}, std::log(2.0));
  near(oracle::gtsam_to_helicoid_vec3(h.transformFrom(Eigen::Vector3d(1, 0, 0))), {1, 4, 3}, "sim3: s R p + x");
}

void sim3_adjoint() {
  // R = R_z(π/2), s = 2, x = (1, 2, 3): helicoid's Ad = [R 0 0; x^R sR −x; 0 0 1], with
  // x^R = [−3 0 2; 0 −3 −1; 1 2 0] and sR = [0 −2 0; 2 0 0; 0 0 2].
  const gtsam::Similarity3 g = oracle::helicoid_to_gtsam_similarity3({kHalf, 0, 0, kHalf}, {1, 2, 3}, std::log(2.0));
  near(oracle::gtsam_to_helicoid_matrix(g.AdjointMap()),
       {0, 1, 0, -3, 0, 1, 0, -1, 0, 0, 0, -3, 2, 0, 0, 0, 1, 2, -1, 0, 0, 0, 0, 0, 0, 2, 0, 0,
        0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, -1, -2, -3, 1},
       "sim3 Ad");
}

void sim3_exp() {
  // [φ; ρ; σ] = [(0, 0, π/2); 0; ln 2] is (q, x, σ) = ((√½, 0, 0, √½), 0, ln 2).
  const auto turn = oracle::gtsam_to_helicoid_sim3(
      gtsam::Similarity3::Expmap(oracle::helicoid_to_gtsam_tangent<7>({0, 0, kPi / 2, 0, 0, 0, std::log(2.0)})));
  near(turn.q, {kHalf, 0, 0, kHalf}, "sim3 exp: q");
  near(turn.x, {0, 0, 0}, "sim3 exp: x");
  near({turn.sigma}, {std::log(2.0)}, "sim3 exp: σ");
  // φ = 0, σ = ln 2: V = (e^σ − 1)/σ I = I/ln 2, so x = ρ/ln 2, which GTSAM holds as t = x/2.
  const double l = std::log(2.0);
  const auto grow = oracle::gtsam_to_helicoid_sim3(
      gtsam::Similarity3::Expmap(oracle::helicoid_to_gtsam_tangent<7>({0, 0, 0, 1, 2, 3, l})));
  near(grow.x, {1 / l, 2 / l, 3 / l}, "sim3 exp: x = V ρ = s t");
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
  rotations();
  tangents();
  matrices();
  se2();
  so3_sides();
  se3();
  extended_pose3();
  sim3_element();
  sim3_adjoint();
  sim3_exp();
  hex_floats();
  if (failures != 0) return 1;
  std::puts("convert_test: all passed");
  return 0;
}
