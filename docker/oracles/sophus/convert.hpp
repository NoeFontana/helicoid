// The corpus's arrays to Sophus's types and back (docs/decisions/0072), at the scalar `T` a record
// is answered in: `double`, or `float` for an `@f32` record (docs/decisions/0077).
//
// helicoid is rotation-first and `w`-first (docs/decisions/0002); Sophus is translation-first:
// SE(2) `[υ; θ]`, SE(3) `[υ; ω]`, Sim(3) `[υ; ω; σ]`. Its Sim(3) stores `√s · q` (`RxSO3`), so the
// scale is read back as `log(‖√s · q‖²)`. Eigen's quaternion constructor takes `(w, x, y, z)` and
// stores `[x, y, z, w]`. The corpus's values are binary64: an input crosses into `T` by `narrow`,
// which refuses to round, and an answer back by the exact widening. Every conversion is a named
// function with a hand-computed test in `convert_test.cpp`, run when the image is built.
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

// The carrier of a record's values and of an answer's: binary64, which holds every binary32.
using Vec = std::vector<double>;

// What Sophus's `SOPHUS_ENSURE` raises (`refused.cpp`): the input is outside its stated domain.
struct Refused : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// --- the scalar ------------------------------------------------------------------------------------

// `x` as a `T`; a binary64 that is not a binary32 is an error, never a rounding.
template <class T>
T narrow(double x) {
  const T y = static_cast<T>(x);
  if (static_cast<double>(y) != x && !std::isnan(x)) {
    throw std::invalid_argument("an input of a binary32 record is not a binary32");
  }
  return y;
}

// --- quaternions -------------------------------------------------------------------------------

template <class T>
Eigen::Quaternion<T> helicoid_to_sophus_quat(const Vec& q) {
  return Eigen::Quaternion<T>(narrow<T>(q.at(0)), narrow<T>(q.at(1)), narrow<T>(q.at(2)), narrow<T>(q.at(3)));
}

template <class T>
Vec sophus_to_helicoid_quat(const Eigen::Quaternion<T>& q) {
  return {static_cast<double>(q.w()), static_cast<double>(q.x()), static_cast<double>(q.y()),
          static_cast<double>(q.z())};
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

template <class T, int N>
Eigen::Matrix<T, N, 1> helicoid_to_sophus_tangent(const Vec& tau, const Order<N>& at) {
  if (tau.size() != static_cast<std::size_t>(N)) {
    throw std::invalid_argument("a tangent has " + std::to_string(N) + " values");
  }
  Eigen::Matrix<T, N, 1> a;
  for (int i = 0; i < N; ++i) a(at[i]) = narrow<T>(tau[i]);
  return a;
}

template <class T, int N>
Vec sophus_to_helicoid_tangent(const Eigen::Matrix<T, N, 1>& a, const Order<N>& at) {
  Vec tau(N);
  for (int i = 0; i < N; ++i) tau[i] = static_cast<double>(a(at[i]));
  return tau;
}

// A Jacobian or an adjoint in Sophus's tangent order as the corpus's array: helicoid's order,
// column-major.
template <class T, int N>
Vec sophus_to_helicoid_matrix(const Eigen::Matrix<T, N, N>& m, const Order<N>& at) {
  Vec out;
  out.reserve(N * N);
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < N; ++i) out.push_back(static_cast<double>(m(at[i], at[j])));
  }
  return out;
}

// A 3 × 3 rotation-block matrix (SO(3)'s Jacobians): the same order, column-major.
template <class T>
Vec sophus_to_helicoid_matrix3(const Eigen::Matrix<T, 3, 3>& m) {
  return Vec(m.data(), m.data() + 9);
}

template <class T>
Eigen::Matrix<T, 3, 1> helicoid_to_sophus_vec3(const Vec& v) {
  return Eigen::Matrix<T, 3, 1>(narrow<T>(v.at(0)), narrow<T>(v.at(1)), narrow<T>(v.at(2)));
}

template <class T>
Vec sophus_to_helicoid_vec3(const Eigen::Matrix<T, 3, 1>& v) {
  return {static_cast<double>(v.x()), static_cast<double>(v.y()), static_cast<double>(v.z())};
}

// --- elements ------------------------------------------------------------------------------------

// helicoid's `z = [cos θ, sin θ]` is Sophus's unit complex `(real, imag)`, normalized by `SO2`.
template <class T>
Sophus::SO2<T> helicoid_to_sophus_so2(const Vec& z) {
  return Sophus::SO2<T>(narrow<T>(z.at(0)), narrow<T>(z.at(1)));
}

template <class T>
Vec sophus_to_helicoid_so2(const Sophus::SO2<T>& r) {
  return {static_cast<double>(r.unit_complex().x()), static_cast<double>(r.unit_complex().y())};
}

template <class T>
Sophus::SE2<T> helicoid_to_sophus_se2(const Vec& z, const Vec& t) {
  return Sophus::SE2<T>(helicoid_to_sophus_so2<T>(z), Eigen::Matrix<T, 2, 1>(narrow<T>(t.at(0)), narrow<T>(t.at(1))));
}

// `SO3`'s quaternion constructor normalizes.
template <class T>
Sophus::SO3<T> helicoid_to_sophus_so3(const Vec& q) {
  return Sophus::SO3<T>(helicoid_to_sophus_quat<T>(q));
}

template <class T>
Sophus::SE3<T> helicoid_to_sophus_se3(const Vec& q, const Vec& x) {
  return Sophus::SE3<T>(helicoid_to_sophus_so3<T>(q), helicoid_to_sophus_vec3<T>(x));
}

// `(q, x, σ)` as `Sim3(RxSO3(e^σ, SO3(q)), x)`: `RxSO3` stores `√(e^σ) · q`, and refuses a scale
// outside `[ε, 1/ε)`, `ε` being Sophus's `Constants<T>::epsilon()`: `10⁻¹⁰` for `double`, `10⁻⁵`
// for `float`. `e^σ` is `T`'s own `exp`.
template <class T>
Sophus::Sim3<T> helicoid_to_sophus_sim3(const Vec& q, const Vec& x, double sigma) {
  return Sophus::Sim3<T>(Sophus::RxSO3<T>(std::exp(narrow<T>(sigma)), helicoid_to_sophus_so3<T>(q)),
                         helicoid_to_sophus_vec3<T>(x));
}

struct HelicoidSim3 {
  Vec q, x;
  double sigma;
};

// `σ = log(scale())`, `scale()` being the stored quaternion's squared norm, as `RxSO3::log` reads it.
template <class T>
HelicoidSim3 sophus_to_helicoid_sim3(const Sophus::Sim3<T>& s) {
  return {sophus_to_helicoid_quat<T>(s.rxso3().so3().unit_quaternion()),
          sophus_to_helicoid_vec3<T>(s.translation()), static_cast<double>(std::log(s.scale()))};
}

}  // namespace oracle
