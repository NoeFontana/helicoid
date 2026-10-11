// The corpus's arrays to Sophus's types and back (docs/decisions/0072).
//
// helicoid is rotation-first and `w`-first (docs/decisions/0002); Sophus is translation-first:
// SE(2) `[υ; θ]`, SE(3) `[υ; ω]`, Sim(3) `[υ; ω; σ]`. Its Sim(3) stores `√s · q` (`RxSO3`), so the
// scale is read back as `log(‖√s · q‖²)`. Eigen's quaternion constructor takes `(w, x, y, z)` and
// stores `[x, y, z, w]`. Every conversion is a named function with a hand-computed test in
// `convert_test.cpp`, run when the image is built.
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <sophus/se2.hpp>
#include <sophus/se3.hpp>
#include <sophus/sim3.hpp>
#include <sophus/so2.hpp>
#include <sophus/so3.hpp>

namespace oracle {

using Vec = std::vector<double>;

// What Sophus's `SOPHUS_ENSURE` raises (`refused.cpp`): the input is outside its stated domain.
struct Refused : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// --- quaternions -------------------------------------------------------------------------------

inline Eigen::Quaterniond helicoid_to_sophus_quat(const Vec& q) {
  return Eigen::Quaterniond(q.at(0), q.at(1), q.at(2), q.at(3));
}

inline Vec sophus_to_helicoid_quat(const Eigen::Quaterniond& q) {
  return {q.w(), q.x(), q.y(), q.z()};
}

// --- tangents and their matrices -----------------------------------------------------------------

// `at[i]`: where helicoid's tangent component `i` sits in Sophus's tangent.
template <int N>
struct Order {
  std::array<int, N> at;
  constexpr int operator[](int i) const { return at[i]; }
};

// SE(2): helicoid `[θ; ρ₀; ρ₁]`, Sophus `[υ₀; υ₁; θ]`.
inline constexpr Order<3> kSe2{{2, 0, 1}};
// SE(3): helicoid `[φ; ρ]`, Sophus `[υ; ω]`.
inline constexpr Order<6> kSe3{{3, 4, 5, 0, 1, 2}};
// Sim(3): helicoid `[φ; ρ; σ]`, Sophus `[υ; ω; σ]`.
inline constexpr Order<7> kSim3{{3, 4, 5, 0, 1, 2, 6}};

template <int N>
Eigen::Matrix<double, N, 1> helicoid_to_sophus_tangent(const Vec& tau, const Order<N>& at) {
  if (tau.size() != static_cast<std::size_t>(N)) {
    throw std::invalid_argument("a tangent has " + std::to_string(N) + " values");
  }
  Eigen::Matrix<double, N, 1> a;
  for (int i = 0; i < N; ++i) a(at[i]) = tau[i];
  return a;
}

template <int N>
Vec sophus_to_helicoid_tangent(const Eigen::Matrix<double, N, 1>& a, const Order<N>& at) {
  Vec tau(N);
  for (int i = 0; i < N; ++i) tau[i] = a(at[i]);
  return tau;
}

// A Jacobian or an adjoint in Sophus's tangent order as the corpus's array: helicoid's order,
// column-major.
template <int N>
Vec sophus_to_helicoid_matrix(const Eigen::Matrix<double, N, N>& m, const Order<N>& at) {
  Vec out;
  out.reserve(N * N);
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < N; ++i) out.push_back(m(at[i], at[j]));
  }
  return out;
}

// A 3 × 3 rotation-block matrix (SO(3)'s Jacobians): the same order, column-major.
inline Vec sophus_to_helicoid_matrix3(const Eigen::Matrix3d& m) {
  return Vec(m.data(), m.data() + 9);
}

inline Eigen::Vector3d helicoid_to_sophus_vec3(const Vec& v) {
  return Eigen::Vector3d(v.at(0), v.at(1), v.at(2));
}

inline Vec sophus_to_helicoid_vec3(const Eigen::Vector3d& v) { return {v.x(), v.y(), v.z()}; }

// --- elements ------------------------------------------------------------------------------------

// helicoid's `z = [cos θ, sin θ]` is Sophus's unit complex `(real, imag)`, normalized by `SO2`.
inline Sophus::SO2d helicoid_to_sophus_so2(const Vec& z) { return Sophus::SO2d(z.at(0), z.at(1)); }

inline Vec sophus_to_helicoid_so2(const Sophus::SO2d& r) {
  return {r.unit_complex().x(), r.unit_complex().y()};
}

inline Sophus::SE2d helicoid_to_sophus_se2(const Vec& z, const Vec& t) {
  return Sophus::SE2d(helicoid_to_sophus_so2(z), Eigen::Vector2d(t.at(0), t.at(1)));
}

// `SO3`'s quaternion constructor normalizes.
inline Sophus::SO3d helicoid_to_sophus_so3(const Vec& q) {
  return Sophus::SO3d(helicoid_to_sophus_quat(q));
}

inline Sophus::SE3d helicoid_to_sophus_se3(const Vec& q, const Vec& x) {
  return Sophus::SE3d(helicoid_to_sophus_so3(q), helicoid_to_sophus_vec3(x));
}

// `(q, x, σ)` as `Sim3(RxSO3(e^σ, SO3(q)), x)`: `RxSO3` stores `√(e^σ) · q`, and refuses a scale
// outside `[10⁻¹⁰, 10¹⁰)`.
inline Sophus::Sim3d helicoid_to_sophus_sim3(const Vec& q, const Vec& x, double sigma) {
  return Sophus::Sim3d(Sophus::RxSO3d(std::exp(sigma), helicoid_to_sophus_so3(q)),
                       helicoid_to_sophus_vec3(x));
}

struct HelicoidSim3 {
  Vec q, x;
  double sigma;
};

// `σ = log(scale())`, `scale()` being the stored quaternion's squared norm, as `RxSO3::log` reads it.
inline HelicoidSim3 sophus_to_helicoid_sim3(const Sophus::Sim3d& s) {
  return {sophus_to_helicoid_quat(s.rxso3().so3().unit_quaternion()),
          sophus_to_helicoid_vec3(s.translation()), std::log(s.scale())};
}

}  // namespace oracle
