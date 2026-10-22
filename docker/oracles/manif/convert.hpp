// The corpus's arrays to manif's types and back (docs/decisions/0075).
//
// helicoid is rotation-first and `w`-first (docs/decisions/0002); manif is translation-first, and
// its SE₂(3) puts the rotation between the two translations: SE(2) `[ρ; θ]`, SE(3) `[ρ; θ]`,
// SE₂(3) `[ρ; θ; ν]` with `ρ` the position's and `ν` the velocity's, where helicoid's SE₂(3) is
// `(R, v, p)` with tangent `[φ; ρ_v; ρ_p]` (docs/maths/se3.md SE.1). Every conversion is a named
// function with a hand-computed test in `convert_test.cpp`, run when the image is built.
#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <manif/manif.h>

// manif's unit-norm checks are its stated domain ("@throws manif::invalid_argument on
// un-normalized quaternion"); the runner answers a record they refuse with `nan`.
#if defined(MANIF_NO_DEBUG) || defined(NDEBUG)
#error "manif's domain checks are the oracle's refusals; build without NDEBUG and MANIF_NO_DEBUG"
#endif

namespace oracle {

using Vec = std::vector<double>;

// manif's world-frame chart on SE(3): `(q Exp φ, x + ρ)`, the product of SO(3) and R³.
using SO3xR3 = manif::Bundle<double, manif::SO3, manif::R3>;

// --- quaternions -------------------------------------------------------------------------------

// Eigen's constructor takes `(w, x, y, z)`; manif stores `[x, y, z, w]`.
inline Eigen::Quaterniond helicoid_to_manif_quat(const Vec& q) {
  return Eigen::Quaterniond(q.at(0), q.at(1), q.at(2), q.at(3));
}

inline Vec manif_to_helicoid_quat(const Eigen::Quaterniond& q) { return {q.w(), q.x(), q.y(), q.z()}; }

// --- tangents and their matrices -----------------------------------------------------------------

// `at[i]`: where helicoid's tangent component `i` sits in manif's tangent.
template <int N>
struct Order {
  std::array<int, N> at;
  constexpr int operator[](int i) const { return at[i]; }
};

// SE(2): helicoid `[θ; ρ₀; ρ₁]`, manif `[ρ₀; ρ₁; θ]`.
inline constexpr Order<3> kSe2{{2, 0, 1}};
// SE(3): helicoid `[φ; ρ]`, manif `[ρ; θ]`.
inline constexpr Order<6> kSe3{{3, 4, 5, 0, 1, 2}};
// SE₂(3): helicoid `[φ; ρ_v; ρ_p]`, manif `[ρ_p; θ; ν]`.
inline constexpr Order<9> kSe23{{3, 4, 5, 6, 7, 8, 0, 1, 2}};

template <int N>
Eigen::Matrix<double, N, 1> helicoid_to_manif_tangent(const Vec& tau, const Order<N>& at) {
  if (tau.size() != static_cast<std::size_t>(N)) {
    throw std::invalid_argument("a tangent has " + std::to_string(N) + " values");
  }
  Eigen::Matrix<double, N, 1> a;
  for (int i = 0; i < N; ++i) a(at[i]) = tau[i];
  return a;
}

template <int N>
Vec manif_to_helicoid_tangent(const Eigen::Matrix<double, N, 1>& a, const Order<N>& at) {
  Vec tau(N);
  for (int i = 0; i < N; ++i) tau[i] = a(at[i]);
  return tau;
}

// A Jacobian or an adjoint in manif's tangent order as the corpus's array: helicoid's order,
// column-major.
template <int N>
Vec manif_to_helicoid_matrix(const Eigen::Matrix<double, N, N>& m, const Order<N>& at) {
  Vec out;
  out.reserve(N * N);
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < N; ++i) out.push_back(m(at[i], at[j]));
  }
  return out;
}

// A 3 × 3 rotation-block matrix (SO(3)'s Jacobians): the same order, column-major.
inline Vec manif_to_helicoid_matrix3(const Eigen::Matrix3d& m) { return Vec(m.data(), m.data() + 9); }

// A `3 x 6` Jacobian of a point by an SE(3) element (0103), manif's `[ρ; θ]` columns as
// helicoid's rotation-first `[φ; ρ]`, column-major.
inline Vec manif_to_helicoid_3x6(const Eigen::Matrix<double, 3, 6>& m) {
  Vec out;
  out.reserve(18);
  for (int j = 0; j < 6; ++j) {
    for (int i = 0; i < 3; ++i) out.push_back(m(i, kSe3[j]));
  }
  return out;
}

inline Eigen::Vector3d helicoid_to_manif_vec3(const Vec& v) { return Eigen::Vector3d(v.at(0), v.at(1), v.at(2)); }

inline Vec manif_to_helicoid_vec3(const Eigen::Vector3d& v) { return {v.x(), v.y(), v.z()}; }

// --- elements ------------------------------------------------------------------------------------

// helicoid's `z = [cos θ, sin θ]` is manif's unit complex `(real, imag)`. manif refuses, never
// normalizes, a non-unit one.
inline manif::SO2d helicoid_to_manif_so2(const Vec& z) { return manif::SO2d(z.at(0), z.at(1)); }

inline Vec manif_to_helicoid_so2(const manif::SO2d& r) { return {r.real(), r.imag()}; }

inline manif::SE2d helicoid_to_manif_se2(const Vec& z, const Vec& t) {
  return manif::SE2d(t.at(0), t.at(1), z.at(0), z.at(1));
}

struct HelicoidSe2 {
  Vec z, t;
};

inline HelicoidSe2 manif_to_helicoid_se2(const manif::SE2d& g) { return {{g.real(), g.imag()}, {g.x(), g.y()}}; }

inline manif::SO3d helicoid_to_manif_so3(const Vec& q) { return manif::SO3d(helicoid_to_manif_quat(q)); }

inline manif::SE3d helicoid_to_manif_se3(const Vec& q, const Vec& x) {
  return manif::SE3d(helicoid_to_manif_vec3(x), helicoid_to_manif_quat(q));
}

// helicoid's `x = [v; p]` is manif's `SE_2_3(t = p, q, v)`.
inline manif::SE_2_3d helicoid_to_manif_se23(const Vec& q, const Vec& x) {
  if (x.size() != 6) throw std::invalid_argument("an SE₂(3) element has 6 translation values");
  return manif::SE_2_3d(Eigen::Vector3d(x[3], x[4], x[5]), helicoid_to_manif_quat(q),
                        Eigen::Vector3d(x[0], x[1], x[2]));
}

struct HelicoidSen3 {
  Vec q, x;
};

inline HelicoidSen3 manif_to_helicoid_se3(const manif::SE3d& g) {
  return {manif_to_helicoid_quat(g.quat()), manif_to_helicoid_vec3(g.translation())};
}

inline HelicoidSen3 manif_to_helicoid_se23(const manif::SE_2_3d& g) {
  const Eigen::Vector3d v = g.linearVelocity(), p = g.translation();
  return {manif_to_helicoid_quat(g.quat()), {v.x(), v.y(), v.z(), p.x(), p.y(), p.z()}};
}

// The world chart's base `(q, x)`; its tangent is `[φ; ρ]`, helicoid's order.
inline SO3xR3 helicoid_to_manif_so3xr3(const Vec& q, const Vec& x) {
  return SO3xR3(helicoid_to_manif_so3(q), manif::R3d(helicoid_to_manif_vec3(x)));
}

inline HelicoidSen3 manif_to_helicoid_so3xr3(const SO3xR3& g) {
  return {manif_to_helicoid_quat(g.element<0>().quat()), manif_to_helicoid_vec3(g.element<1>().coeffs())};
}

}  // namespace oracle
