// Hand-computed tests of every conversion in `convert.hpp` and of the hex floats, on values whose
// components all differ, so a permuted, transposed or reversed layout fails here rather than scoring
// (docs/PHASE1.md §7), at `double` and at `float` (docs/decisions/0077). The image build runs this
// and stops on a failure.
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

// Far above the rounding of the inputs and of Sophus at these values, far below any slip: `double`'s,
// then `float`'s.
constexpr double kTol = 1e-14;
constexpr double kTolF = 1e-6;

void near(const oracle::Vec& got, const oracle::Vec& want, const std::string& what, double tol = kTol) {
  bool ok = got.size() == want.size();
  for (std::size_t i = 0; ok && i < got.size(); ++i) ok = std::fabs(got[i] - want[i]) <= tol;
  check(ok, what);
}

void exact(const oracle::Vec& got, const oracle::Vec& want, const std::string& what) {
  check(got == want, what);
}

bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

// Every value is a binary32: the widening of a `float`, never a `double` result.
bool binary32(const oracle::Vec& v) {
  for (const double x : v) {
    if (!std::isnan(x) && static_cast<double>(static_cast<float>(x)) != x) return false;
  }
  return true;
}

// `x`'s nearest binary32, as the corpus would hold it in an `@f32` record.
double f32(double x) { return static_cast<double>(static_cast<float>(x)); }

void quaternions() {
  const Eigen::Quaterniond q = oracle::helicoid_to_sophus_quat<double>({0.1, 0.2, 0.3, 0.4});
  check(q.w() == 0.1 && q.x() == 0.2 && q.y() == 0.3 && q.z() == 0.4, "quat: (w, x, y, z)");
  check(q.coeffs()(3) == 0.1 && q.coeffs()(0) == 0.2, "quat: Eigen stores w last");
  exact(oracle::sophus_to_helicoid_quat(q), {0.1, 0.2, 0.3, 0.4}, "quat: back");
  const Eigen::Quaternionf p = oracle::helicoid_to_sophus_quat<float>({0.5, 0.25, -0.125, 0.75});
  check(p.w() == 0.5f && p.x() == 0.25f && p.y() == -0.125f && p.z() == 0.75f, "quat f32: (w, x, y, z)");
  check(p.coeffs()(3) == 0.5f && p.coeffs()(0) == 0.25f, "quat f32: Eigen stores w last");
  exact(oracle::sophus_to_helicoid_quat(p), {0.5, 0.25, -0.125, 0.75}, "quat f32: back");
}

void narrowing() {
  // 0.1 is a binary64 and no binary32; its nearest binary32 is 13421773 · 2⁻²⁷ = 0x1.99999ap-4.
  bool refused = false;
  try {
    oracle::narrow<float>(0.1);
  } catch (const std::invalid_argument&) {
    refused = true;
  }
  check(refused, "narrow: 0.1 is not a binary32, refused");
  check(oracle::narrow<float>(f32(0.1)) == 0.1f, "narrow: 0x1.99999ap-4 is");
  check(oracle::narrow<double>(0.1) == 0.1, "narrow: anything is a binary64");
  check(oracle::narrow<float>(1.0 + std::ldexp(1.0, -23)) == 1.0f + 0x1p-23f, "narrow: one ulp above 1");
  check(std::isnan(oracle::narrow<float>(std::nan(""))), "narrow: nan");
  check(oracle::format_hex(f32(0.1)) == "0x1.99999ap-4", "hex: a binary32 widened prints its own digits");
}

void tangents() {
  using oracle::helicoid_to_sophus_tangent;
  using oracle::sophus_to_helicoid_tangent;
  const auto se2 = helicoid_to_sophus_tangent<double>({1, 2, 3}, oracle::kSe2);
  check(se2 == Eigen::Vector3d(2, 3, 1), "se2 tangent: [θ; ρ] → [υ; θ]");
  exact(sophus_to_helicoid_tangent(se2, oracle::kSe2), {1, 2, 3}, "se2 tangent: back");
  const auto se3 = helicoid_to_sophus_tangent<double>({1, 2, 3, 4, 5, 6}, oracle::kSe3);
  check(se3 == (Sophus::SE3d::Tangent() << 4, 5, 6, 1, 2, 3).finished(), "se3 tangent: [φ; ρ] → [υ; ω]");
  exact(sophus_to_helicoid_tangent(se3, oracle::kSe3), {1, 2, 3, 4, 5, 6}, "se3 tangent: back");
  const auto sim3 = helicoid_to_sophus_tangent<double>({1, 2, 3, 4, 5, 6, 7}, oracle::kSim3);
  check(sim3 == (Sophus::Sim3d::Tangent() << 4, 5, 6, 1, 2, 3, 7).finished(),
        "sim3 tangent: [φ; ρ; σ] → [υ; ω; σ]");
  exact(sophus_to_helicoid_tangent(sim3, oracle::kSim3), {1, 2, 3, 4, 5, 6, 7}, "sim3 tangent: back");
  const auto sim3f = helicoid_to_sophus_tangent<float>({1, 2, 3, 4, 5, 6, 7}, oracle::kSim3);
  check(sim3f == (Sophus::Sim3f::Tangent() << 4, 5, 6, 1, 2, 3, 7).finished(),
        "sim3 tangent f32: [φ; ρ; σ] → [υ; ω; σ]");
  exact(sophus_to_helicoid_tangent(sim3f, oracle::kSim3), {1, 2, 3, 4, 5, 6, 7}, "sim3 tangent f32: back");
}

void matrices() {
  // m(r, c) = 10 r + c in Sophus's order [υ₀; υ₁; θ]; helicoid's (i, j) is Sophus's (at[i], at[j]).
  Eigen::Matrix3d m;
  m << 0, 1, 2, 10, 11, 12, 20, 21, 22;
  exact(oracle::sophus_to_helicoid_matrix(m, oracle::kSe2), {22, 2, 12, 20, 0, 10, 21, 1, 11},
        "se2 matrix: permuted, column-major");
  exact(oracle::sophus_to_helicoid_matrix3(m), {0, 10, 20, 1, 11, 21, 2, 12, 22}, "matrix3: column-major");
  const Eigen::Matrix3f mf = m.cast<float>();
  exact(oracle::sophus_to_helicoid_matrix(mf, oracle::kSe2), {22, 2, 12, 20, 0, 10, 21, 1, 11},
        "se2 matrix f32: permuted, column-major");
}

void se2_adjoint() {
  // θ = π/2, t = (1, 2): helicoid's Ad = [1 0 0; t₁ c −s; −t₀ s c] = [1 0 0; 2 0 −1; −1 1 0].
  const Sophus::SE2d g = oracle::helicoid_to_sophus_se2<double>({0, 1}, {1, 2});
  near(oracle::sophus_to_helicoid_matrix(g.Adj(), oracle::kSe2), {1, 2, -1, 0, 0, 1, 0, -1, 0}, "se2 Ad");
  near(oracle::sophus_to_helicoid_so2(oracle::helicoid_to_sophus_so2<double>({0.6, 0.8})), {0.6, 0.8},
       "so2: (re, im)");
}

void se3_left_jacobian() {
  // φ = 0, ρ = (1, 2, 3): helicoid's J_l = [I 0; ½ρ^ I], ½ρ^ = [0 −1.5 1; 1.5 0 −0.5; −1 0.5 0].
  const auto tau = oracle::helicoid_to_sophus_tangent<double>({0, 0, 0, 1, 2, 3}, oracle::kSe3);
  exact(oracle::sophus_to_helicoid_matrix(Sophus::SE3d::leftJacobian(tau), oracle::kSe3),
        {1, 0, 0, 0, 1.5, -1, 0, 1, 0, -1.5, 0, 0.5, 0, 0, 1, 1, -0.5, 0,
         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1},
        "se3 J_l at φ = 0");
}

void se3_element() {
  const double h = std::sqrt(0.5);
  const Sophus::SE3d g = oracle::helicoid_to_sophus_se3<double>({h, 0, 0, h}, {1, 2, 3});
  near(oracle::sophus_to_helicoid_vec3<double>(g * Eigen::Vector3d(1, 0, 0)), {1, 3, 3}, "se3: R_z(π/2) e_x + t");
}

void sim3_element() {
  // q = (1, 2, 3, 4)/√30, σ = ln 4: Sophus stores 2q, and reads σ back from its squared norm.
  const double r = std::sqrt(30.0);
  const oracle::Vec q{1 / r, 2 / r, 3 / r, 4 / r};
  const Sophus::Sim3d x = oracle::helicoid_to_sophus_sim3<double>(q, {5, 6, 7}, std::log(4.0));
  const Eigen::Quaterniond& stored = x.rxso3().quaternion();
  near({stored.w(), stored.x(), stored.y(), stored.z()}, {2 / r, 4 / r, 6 / r, 8 / r}, "sim3: √s · q");
  near(oracle::sophus_to_helicoid_vec3<double>(x.translation()), {5, 6, 7}, "sim3: x");
  const oracle::HelicoidSim3 back = oracle::sophus_to_helicoid_sim3(x);
  near(back.q, q, "sim3: q back");
  near({back.sigma}, {std::log(4.0)}, "sim3: σ back");
  bool refused = false;
  try {
    oracle::helicoid_to_sophus_sim3<double>(q, {5, 6, 7}, 30.0);
  } catch (const oracle::Refused&) {
    refused = true;
  }
  check(refused, "sim3: e^30 is past RxSO3's 10^10, refused");
}

void sim3_element_f32() {
  // q = (1, 1, 1, 1)/2 and σ = fl₃₂(ln 4), each exactly a binary32: Sophus stores √(e^σ) q, which
  // is 2q = (1, 1, 1, 1) to float's rounding, and reads σ back to it; every value it gives is a
  // binary32.
  const oracle::Vec q{0.5, 0.5, 0.5, 0.5};
  const double sigma = f32(std::log(4.0));
  const Sophus::Sim3f x = oracle::helicoid_to_sophus_sim3<float>(q, {5, 6, 7}, sigma);
  const Eigen::Quaternionf& stored = x.rxso3().quaternion();
  const oracle::Vec two_q{stored.w(), stored.x(), stored.y(), stored.z()};
  near(two_q, {1, 1, 1, 1}, "sim3 f32: √s · q", kTolF);
  check(binary32(two_q), "sim3 f32: the stored quaternion is binary32");
  const oracle::HelicoidSim3 back = oracle::sophus_to_helicoid_sim3(x);
  near(back.q, q, "sim3 f32: q back", kTolF);
  exact(back.x, {5, 6, 7}, "sim3 f32: x back");
  near({back.sigma}, {std::log(4.0)}, "sim3 f32: σ back", kTolF);
  check(binary32(back.q) && binary32({back.sigma}), "sim3 f32: back is binary32");
  // `float`'s domain is Sophus's `ε = 10⁻⁵`: e^12 ≈ 1.6 · 10⁵ is refused at `float` and not at
  // `double`.
  bool refused = false;
  try {
    oracle::helicoid_to_sophus_sim3<float>(q, {5, 6, 7}, 12.0);
  } catch (const oracle::Refused&) {
    refused = true;
  }
  check(refused, "sim3 f32: e^12 is past RxSO3<float>'s 10^5, refused");
  oracle::helicoid_to_sophus_sim3<double>(q, {5, 6, 7}, 12.0);
  // An input that is not a binary32 is an error, not a rounding.
  bool rounded = false;
  try {
    oracle::helicoid_to_sophus_sim3<float>(q, {0.1, 6, 7}, sigma);
  } catch (const std::invalid_argument&) {
    rounded = true;
  }
  check(rounded, "sim3 f32: x = 0.1 is not a binary32, refused");
}

void sim3_adjoint() {
  // R = R_z(π/2), s = 2, t = (1, 2, 3): helicoid's Ad = [R 0 0; t^R sR −t; 0 0 1], with
  // t^R = [−3 0 2; 0 −3 −1; 1 2 0] and sR = [0 −2 0; 2 0 0; 0 0 2].
  const oracle::Vec want{0, 1, 0, -3, 0, 1, 0, -1, 0, 0, 0, -3, 2, 0, 0, 0, 1, 2, -1, 0, 0, 0, 0, 0, 0, 2, 0, 0,
                         0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, -1, -2, -3, 1};
  const double h = std::sqrt(0.5);
  const Sophus::Sim3d x = oracle::helicoid_to_sophus_sim3<double>({h, 0, 0, h}, {1, 2, 3}, std::log(2.0));
  near(oracle::sophus_to_helicoid_matrix(x.Adj(), oracle::kSim3), want, "sim3 Ad");
  // The same at `float`, from the binary32 nearest √½ and ln 2.
  const Sophus::Sim3f xf =
      oracle::helicoid_to_sophus_sim3<float>({f32(h), 0, 0, f32(h)}, {1, 2, 3}, f32(std::log(2.0)));
  const oracle::Vec ad = oracle::sophus_to_helicoid_matrix(xf.Adj(), oracle::kSim3);
  near(ad, want, "sim3 Ad f32", kTolF);
  check(binary32(ad), "sim3 Ad f32: binary32");
}

void sim3_exp() {
  // [φ; ρ; σ] = [(0, 0, π/2); 0; ln 2] is (q, x, σ) = ((√½, 0, 0, √½), 0, ln 2); with φ = σ = 0,
  // x = ρ.
  const double h = std::sqrt(0.5);
  const auto turn = oracle::sophus_to_helicoid_sim3(Sophus::Sim3d::exp(oracle::helicoid_to_sophus_tangent<double>(
      {0, 0, 2 * std::atan(1.0), 0, 0, 0, std::log(2.0)}, oracle::kSim3)));
  near(turn.q, {h, 0, 0, h}, "sim3 exp: q");
  near(turn.x, {0, 0, 0}, "sim3 exp: x");
  near({turn.sigma}, {std::log(2.0)}, "sim3 exp: σ");
  const auto shift = oracle::sophus_to_helicoid_sim3(
      Sophus::Sim3d::exp(oracle::helicoid_to_sophus_tangent<double>({0, 0, 0, 1, 2, 3, 0}, oracle::kSim3)));
  near(shift.x, {1, 2, 3}, "sim3 exp: x = ρ");
  // At `float`: the same turn from binary32 inputs, every output a binary32; and the shift exact.
  const auto turn_f = oracle::sophus_to_helicoid_sim3(Sophus::Sim3f::exp(oracle::helicoid_to_sophus_tangent<float>(
      {0, 0, f32(2 * std::atan(1.0)), 0, 0, 0, f32(std::log(2.0))}, oracle::kSim3)));
  near(turn_f.q, {h, 0, 0, h}, "sim3 exp f32: q", kTolF);
  near(turn_f.x, {0, 0, 0}, "sim3 exp f32: x", kTolF);
  near({turn_f.sigma}, {std::log(2.0)}, "sim3 exp f32: σ", kTolF);
  check(binary32(turn_f.q) && binary32(turn_f.x) && binary32({turn_f.sigma}), "sim3 exp f32: binary32");
  const auto shift_f = oracle::sophus_to_helicoid_sim3(
      Sophus::Sim3f::exp(oracle::helicoid_to_sophus_tangent<float>({0, 0, 0, 1, 2, 3, 0}, oracle::kSim3)));
  exact(shift_f.x, {1, 2, 3}, "sim3 exp f32: x = ρ");
}

// Column `i` of `got` against `want`'s, up to its sign: an eigenvector's sign is the solver's.
void near_columns(const oracle::Vec& got, const oracle::Vec& want, const std::string& what, double tol) {
  bool ok = got.size() == 9 && want.size() == 9;
  for (std::size_t c = 0; ok && c < 3; ++c) {
    bool same = true, flipped = true;
    for (std::size_t r = 0; r < 3; ++r) {
      same = same && std::fabs(got[3 * c + r] - want[3 * c + r]) <= tol;
      flipped = flipped && std::fabs(got[3 * c + r] + want[3 * c + r]) <= tol;
    }
    ok = same || flipped;
  }
  check(ok, what);
}

void eig3() {
  // The corpus's `A` is column-major, as Eigen's storage: entry k is (k mod 3, k div 3).
  const Eigen::Matrix3d m = oracle::helicoid_to_eigen_matrix3<double>({0, 1, 2, 3, 4, 5, 6, 7, 8});
  check(m(1, 0) == 1 && m(0, 1) == 3 && m(2, 1) == 5, "eig3: A column-major");
  // A = [5 0 0; 0 2 1; 0 1 2], its largest eigenvalue first on the diagonal: λ = (1, 3, 5)
  // ascending, V's columns (0, 1, −1)/√2, (0, 1, 1)/√2, e₀; a transposed V would fail.
  const oracle::Vec a{5, 0, 0, 0, 2, 1, 0, 1, 2};
  const double h = std::sqrt(0.5);
  const oracle::Vec v{0, h, -h, 0, h, h, 1, 0, 0};
  using Solver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>;
  const oracle::HelicoidEig3 e = oracle::eigen_to_helicoid_eig3<double>(Solver(oracle::helicoid_to_eigen_matrix3<double>(a)));
  near(e.lambda, {1, 3, 5}, "eig3: λ ascending");
  near_columns(e.v, v, "eig3: V's columns, column-major", kTol);
  using SolverF = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f>;
  const oracle::HelicoidEig3 f = oracle::eigen_to_helicoid_eig3<float>(SolverF(oracle::helicoid_to_eigen_matrix3<float>(a)));
  near(f.lambda, {1, 3, 5}, "eig3 f32: λ ascending", kTolF);
  near_columns(f.v, v, "eig3 f32: V's columns, column-major", kTolF);
  check(binary32(f.lambda) && binary32(f.v), "eig3 f32: binary32");
}

// `U diag(σ) Vᵀ` from `svd3`'s column-major arrays, against `a` (column-major), and `U, V ∈ SO(3)`
// and `σ₁ ≥ σ₂ ≥ |σ₃|`, each to `tol`.
void signed_svd(const oracle::HelicoidSvd3& s, const oracle::Vec& a, const std::string& what, double tol) {
  const Eigen::Map<const Eigen::Matrix3d> u(s.u.data()), v(s.v.data()), want(a.data());
  const Eigen::Vector3d sigma(s.sigma[0], s.sigma[1], s.sigma[2]);
  check((u * sigma.asDiagonal() * v.transpose() - want).norm() <= tol, what + ": U Σ Vᵀ = A");
  check(std::fabs(u.determinant() - 1) <= tol && std::fabs(v.determinant() - 1) <= tol, what + ": det U = det V = 1");
  check((u.transpose() * u - Eigen::Matrix3d::Identity()).norm() <= tol &&
            (v.transpose() * v - Eigen::Matrix3d::Identity()).norm() <= tol,
        what + ": U, V orthogonal");
  check(sigma(0) >= sigma(1) && sigma(1) >= std::fabs(sigma(2)), what + ": σ₁ ≥ σ₂ ≥ |σ₃|");
}

template <class T>
oracle::HelicoidSvd3 svd3_of(const oracle::Vec& a) {
  return oracle::eigen_to_helicoid_svd3<T>(Eigen::JacobiSVD<Eigen::Matrix<T, 3, 3>>(
      oracle::helicoid_to_eigen_matrix3<T>(a), Eigen::ComputeFullU | Eigen::ComputeFullV));
}

void svd3() {
  // diag(2, −3, 1), det < 0: σ = (3, 2, −1), the negative one the smallest in magnitude.
  const oracle::Vec mixed{2, 0, 0, 0, -3, 0, 0, 0, 1};
  const oracle::HelicoidSvd3 d = svd3_of<double>(mixed);
  near(d.sigma, {3, 2, -1}, "svd3: diag(2, −3, 1)'s σ");
  signed_svd(d, mixed, "svd3: diag(2, −3, 1)", kTol);
  // diag(1, 1, −5): σ = (5, 1, −1), ordered by magnitude before it is signed (0015 PH.2).
  const oracle::Vec order{1, 0, 0, 0, 1, 0, 0, 0, -5};
  const oracle::HelicoidSvd3 o = svd3_of<double>(order);
  near(o.sigma, {5, 1, -1}, "svd3: diag(1, 1, −5)'s σ");
  signed_svd(o, order, "svd3: diag(1, 1, −5)", kTol);
  // R_z(π/2): σ = (1, 1, 1) and U Vᵀ = R.
  const oracle::Vec turn{0, 1, 0, -1, 0, 0, 0, 0, 1};
  const oracle::HelicoidSvd3 r = svd3_of<double>(turn);
  near(r.sigma, {1, 1, 1}, "svd3: a rotation's σ");
  const Eigen::Map<const Eigen::Matrix3d> ru(r.u.data()), rv(r.v.data()), rr(turn.data());
  check((ru * rv.transpose() - rr).norm() <= kTol, "svd3: a rotation's U Vᵀ");
  // [[1, 2, 3], [4, 5, 6], [7, 8, 10]], det = −3, column-major: not symmetric, so a transposed U
  // or V fails the reconstruction; and its negation, det = 3.
  const oracle::Vec general{1, 4, 7, 2, 5, 8, 3, 6, 10};
  const oracle::HelicoidSvd3 g = svd3_of<double>(general);
  signed_svd(g, general, "svd3: general, det < 0", 1e-13);
  check(g.sigma[2] < 0, "svd3: sign σ₃ = sign det A < 0");
  oracle::Vec negated(general);
  for (double& x : negated) x = -x;
  const oracle::HelicoidSvd3 n = svd3_of<double>(negated);
  signed_svd(n, negated, "svd3: general, det > 0", 1e-13);
  check(n.sigma[2] > 0, "svd3: sign σ₃ = sign det A > 0");
  // At `float`: the same, every output a binary32.
  const oracle::HelicoidSvd3 f = svd3_of<float>(general);
  signed_svd(f, general, "svd3 f32: general", 1e-5);
  check(f.sigma[2] < 0, "svd3 f32: sign σ₃");
  check(binary32(f.u) && binary32(f.sigma) && binary32(f.v), "svd3 f32: binary32");
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
  narrowing();
  tangents();
  matrices();
  se2_adjoint();
  se3_left_jacobian();
  se3_element();
  sim3_element();
  sim3_element_f32();
  sim3_adjoint();
  sim3_exp();
  eig3();
  svd3();
  hex_floats();
  if (failures != 0) return 1;
  std::puts("convert_test: all passed");
  return 0;
}
