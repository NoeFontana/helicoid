// Hand-computed tests of every conversion in `convert.hpp` and of the hex floats, on values whose
// components all differ, so a permuted, transposed or reversed layout fails here rather than scoring
// (docs/PHASE1.md §7); then each id of `spline.hpp` against what it must equal: the Jacobians
// against central differences of the helper's own value and rates under helicoid's right
// perturbation `X_m Exp(δ)`, which settles their side and sign (docs/NUMERICS.md §21), the twists'
// order on a translation-only window, and the Gram matrix against §21's check values. The image
// build runs this and stops on a failure.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

#include "convert.hpp"
#include "hexfloat.hpp"
#include "spline.hpp"

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", what.c_str());
    ++failures;
  }
}

constexpr double kTol = 1e-14;

void near(const oracle::Vec& got, const oracle::Vec& want, const std::string& what, double tol = kTol) {
  bool ok = got.size() == want.size();
  for (std::size_t i = 0; ok && i < got.size(); ++i) ok = std::fabs(got[i] - want[i]) <= tol;
  check(ok, what);
}

void exact(const oracle::Vec& got, const oracle::Vec& want, const std::string& what) { check(got == want, what); }

bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

bool refused_off_grid(double u, double inv_h) {
  try {
    oracle::nanoseconds(u, inv_h);
  } catch (const oracle::OffGrid&) {
    return true;
  }
  return false;
}

void elements() {
  const oracle::Vec q = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
  const auto r = oracle::sophus_so3_params(q, 1);
  check(r == std::array<double, 4>{0.6, 0.7, 0.8, 0.5}, "so3: control 1 as [x, y, z, w]");
  const Sophus::SO3d mapped = oracle::sophus_so3(r);
  check(mapped.unit_quaternion().w() == 0.5 && mapped.unit_quaternion().x() == 0.6, "so3: mapped, not normalized");
  exact(oracle::sophus_to_helicoid_quat(mapped.unit_quaternion()), {0.5, 0.6, 0.7, 0.8}, "quat: back w-first");
  const auto x = oracle::sophus_se3_params(q, {1, 2, 3, 4, 5, 6}, 1);
  check(x == std::array<double, 7>{0.6, 0.7, 0.8, 0.5, 4, 5, 6}, "se3: control 1 as quaternion, translation");
  const Sophus::SE3d pose = Eigen::Map<const Sophus::SE3d>(x.data());
  check(pose.translation() == Eigen::Vector3d(4, 5, 6) && pose.unit_quaternion().w() == 0.5,
        "se3: Sophus's storage reads as written");
  Sophus::SE3d::Tangent a;
  a << 1, 2, 3, 4, 5, 6;
  exact(oracle::sophus_to_helicoid_twist(a), {4, 5, 6, 1, 2, 3}, "twist: [υ; ω] to [ω; v]");
  Eigen::Matrix3d m;
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  exact(oracle::column_major(m), {1, 4, 7, 2, 5, 8, 3, 6, 9}, "matrix: column-major");
}

void sides() {
  // Quarter turns: `R_z(π/2)` maps x to y, `R_x(π/2)` maps y to z, exactly.
  Eigen::Matrix3d rz, rx, left;
  rz << 0, -1, 0, 1, 0, 0, 0, 0, 1;
  rx << 1, 0, 0, 0, 0, -1, 0, 1, 0;
  left << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  const Sophus::SO3d value(rz), control(rx);
  // `L R_x`: columns `L e₀, L e₂, −L e₁`.
  Eigen::Matrix3d rate;
  rate << 1, 3, -2, 4, 6, -5, 7, 9, -8;
  near(oracle::column_major(oracle::right_rate_jacobian(left, control)), oracle::column_major(rate), "J_ω: J^L R_m");
  // `R_zᵀ (L R_x)`: rows `y, −x, z` of `L R_x`.
  Eigen::Matrix3d val;
  val << 4, 6, -5, -1, -3, 2, 7, 9, -8;
  near(oracle::column_major(oracle::right_value_jacobian(left, value, control)), oracle::column_major(val),
       "J_X: R_Xᵀ J^L R_m");
}

void time_grid() {
  const auto ns = [](double u, double inv_h) {
    const auto t = oracle::nanoseconds(u, inv_h);
    return std::array<long long, 2>{t.dt, t.t};
  };
  check(ns(0, 1) == std::array<long long, 2>{1000000000, 0}, "ns: u = 0, h = 1 s");
  check(ns(0.5, 1) == std::array<long long, 2>{1000000000, 500000000}, "ns: u = ½");
  check(ns(0.5, 1000) == std::array<long long, 2>{1000000, 500000}, "ns: h = 1 ms");
  check(ns(0.25, 0.001) == std::array<long long, 2>{1000000000000, 250000000000}, "ns: h = 1000 s");
  check(ns(0, 0x1p-10) == std::array<long long, 2>{1024000000000, 0}, "ns: h = 2¹⁰ s");
  check(refused_off_grid(0.1, 1), "ns: 0.1 · 10⁹ is no integer");
  check(refused_off_grid(1 - 0x1p-53, 1), "ns: 1 − 2⁻⁵³ is no integer");
  check(refused_off_grid(0, 1024), "ns: h = 2⁻¹⁰ s is 976562.5 ns");
  check(refused_off_grid(0.5, 3), "ns: 10⁹ / 3 is no integer");
}

void scaling() {
  exact(oracle::scaled({2, 4}, 0, 8), {2, 4}, "scaled: a value is not");
  exact(oracle::scaled({2, 4}, 1, 8), {0.25, 0.5}, "scaled: ω / h⁻¹");
  exact(oracle::scaled({2, 4}, 2, 8), {2.0 / 64, 4.0 / 64}, "scaled: ω̇ / h⁻²");
}

// An SO(3) window of order 4 with increments well inside `Log`'s domain, at `u = ½`, `h = 1`.
oracle::Fields so3_window() {
  oracle::Fields in{{"u", {0.5}}, {"inv_h", {1}}};
  oracle::Vec& q = in["q"];
  const double phi[4][3] = {{0.1, -0.2, 0.3}, {0.4, 0.1, -0.5}, {-0.3, 0.7, 0.2}, {0.6, -0.4, 0.9}};
  for (const auto& p : phi) {
    const auto v = oracle::sophus_to_helicoid_quat(Sophus::SO3d::exp(Eigen::Vector3d(p[0], p[1], p[2])).unit_quaternion());
    q.insert(q.end(), v.begin(), v.end());
  }
  return in;
}

// `in` with control `m` moved to `X_m Exp(e)`.
oracle::Fields perturbed(oracle::Fields in, int m, const Eigen::Vector3d& e) {
  oracle::Vec& q = in["q"];
  const Sophus::SO3d x = oracle::sophus_so3(oracle::sophus_so3_params(q, m)) * Sophus::SO3d::exp(e);
  const auto v = oracle::sophus_to_helicoid_quat(x.unit_quaternion());
  std::copy(v.begin(), v.end(), q.begin() + 4 * m);
  return in;
}

// The output of order `D` as a tangent: the value's `Log(X⁻¹ X')` against `base`, a rate itself.
template <int D>
Eigen::Vector3d so3_output(const oracle::Fields& in, const Sophus::SO3d& base) {
  const oracle::Fields out = oracle::so3_rate<D>(in);
  if constexpr (D == 0) {
    const oracle::Vec q = out.at("q");
    const Sophus::SO3d x(Eigen::Quaterniond(q[0], q[1], q[2], q[3]));
    return (base.inverse() * x).log();
  } else {
    const oracle::Vec v = out.at(D == 1 ? "omega" : "alpha");
    return Eigen::Vector3d(v[0], v[1], v[2]);
  }
}

template <int D>
void so3_jacobian_is_right() {
  const oracle::Fields in = so3_window();
  const oracle::Vec q0 = oracle::so3_rate<0>(in).at("q");
  const Sophus::SO3d base(Eigen::Quaterniond(q0[0], q0[1], q0[2], q0[3]));
  const oracle::Vec jac = oracle::so3_jacobian<D>(in).at("jac");
  const double eps = 1e-6;
  oracle::Vec fd;
  for (int m = 0; m < 4; ++m) {
    for (int k = 0; k < 3; ++k) {
      const Eigen::Vector3d e = Eigen::Vector3d::Unit(k) * eps;
      const Eigen::Vector3d d =
          (so3_output<D>(perturbed(in, m, e), base) - so3_output<D>(perturbed(in, m, -e), base)) / (2 * eps);
      fd.insert(fd.end(), {d.x(), d.y(), d.z()});
    }
  }
  near(jac, fd, "so3 Jacobian of order " + std::to_string(D) + ": right, against central differences", 1e-8);
}

void rn_jacobian_is_the_weight() {
  oracle::Fields in{{"u", {0.5}}, {"inv_h", {1}}, {"p", {}}};
  for (int i = 0; i < 5; ++i) in["p"].insert(in["p"].end(), {1.0 * i, 0.5 * i * i, -2.0 * i});
  for (const auto& [d, field] : {std::pair{0, "p"}, {1, "v"}, {2, "a"}}) {
    const oracle::Vec k = d == 0 ? oracle::rn_jacobian<0>(in).at("k")
                          : d == 1 ? oracle::rn_jacobian<1>(in).at("k")
                                   : oracle::rn_jacobian<2>(in).at("k");
    for (int m = 0; m < 5; ++m) {
      oracle::Fields unit = in;
      std::fill(unit["p"].begin(), unit["p"].end(), 0.0);
      unit["p"][3 * m] = 1;
      const oracle::Fields out = d == 0 ? oracle::rn_rate<0>(unit) : d == 1 ? oracle::rn_rate<1>(unit) : oracle::rn_rate<2>(unit);
      near({k[m]}, {out.at(field)[0]}, "rn Jacobian of order " + std::to_string(d) + ": control " + std::to_string(m), 1e-14);
    }
  }
}

// Collinear, evenly spaced, unrotated controls: the twist is `[0; Δ]` per unit time, its
// derivative 0, and Rⁿ's velocity `Δ`.
void twists_are_rotation_first() {
  const Eigen::Vector3d step(1, 2, 3);
  oracle::Fields se3{{"u", {0.25}}, {"inv_h", {1}}, {"q", {}}, {"x", {}}};
  oracle::Fields rn{{"u", {0.25}}, {"inv_h", {1}}, {"p", {}}};
  for (int i = 0; i < 4; ++i) {
    se3["q"].insert(se3["q"].end(), {1, 0, 0, 0});
    for (int c = 0; c < 3; ++c) {
      se3["x"].push_back(i * step(c));
      rn["p"].push_back(i * step(c));
    }
  }
  near(oracle::se3_rate<1>(se3).at("twist"), {0, 0, 0, 1, 2, 3}, "se3 twist: [ω; v]", 1e-13);
  near(oracle::se3_rate<2>(se3).at("twist"), {0, 0, 0, 0, 0, 0}, "se3 twist derivative", 1e-13);
  near(oracle::rn_rate<1>(rn).at("v"), {1, 2, 3}, "rn velocity", 1e-13);
  near(oracle::se3_rate<0>(se3).at("x"), oracle::rn_rate<0>(rn).at("p"), "se3 value: x as Rⁿ's p", 1e-13);
}

// docs/NUMERICS.md §21's check values at `h = 1`: `G⁽¹⁾₀₀ = 1/20` at `K = 4`, `G⁽²⁾₀₀ = 1/252` at
// `K = 6`; every row sums to 0.
void gram() {
  const oracle::Vec g4 = oracle::rn_gram({{"k", {4}}, {"r", {1}}, {"inv_h", {1}}}).at("g");
  const oracle::Vec g6 = oracle::rn_gram({{"k", {6}}, {"r", {2}}, {"inv_h", {1}}}).at("g");
  near({g4[0]}, {1.0 / 20}, "gram: G⁽¹⁾₀₀ at K = 4");
  near({g6[0]}, {1.0 / 252}, "gram: G⁽²⁾₀₀ at K = 6");
  for (int s = 0; s < 6; ++s) {
    double sum = 0;
    for (int t = 0; t < 6; ++t) sum += g6[6 * t + s];
    near({sum}, {0}, "gram: row " + std::to_string(s) + " sums to 0");
  }
  const oracle::Vec slow = oracle::rn_gram({{"k", {4}}, {"r", {1}}, {"inv_h", {0x1p-10}}}).at("g");
  near({slow[0]}, {0x1p-10 / 20}, "gram: G⁽¹⁾ scales by h⁻¹");
}

void hex_floats() {
  check(oracle::format_hex(1.0) == "0x1p+0", "hex: 1");
  check(oracle::format_hex(-0.0) == "-0x0p+0", "hex: −0");
  check(oracle::format_hex(std::nan("")) == "nan" && oracle::format_hex(-INFINITY) == "-inf", "hex: words");
  check(same_bits(oracle::parse_hex("-0x0.0p+0"), -0.0), "hex: −0 in");
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
  elements();
  sides();
  time_grid();
  scaling();
  so3_jacobian_is_right<0>();
  so3_jacobian_is_right<1>();
  so3_jacobian_is_right<2>();
  rn_jacobian_is_the_weight();
  twists_are_rotation_first();
  gram();
  hex_floats();
  if (failures != 0) return 1;
  std::puts("convert_test: all passed");
  return 0;
}
