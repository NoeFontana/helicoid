// The corpus's spline records to basalt-headers' types and back (docs/decisions/0116 decision 4).
//
// basalt's knots are Sophus elements: SO(3) in Eigen's quaternion storage `[x, y, z, w]`, SE(3)
// that and then the translation, its tangent `[υ; ω]`; helicoid is `w`-first and rotation-first
// (docs/decisions/0002). A knot is read through `Eigen::Map`, as `CeresSplineHelper` reads its
// parameter blocks, so neither path normalizes what the record holds. `So3Spline`'s Jacobians are
// left on both sides (its own tests perturb `Exp(x) X_m` and read `Log(X' X⁻¹)`); helicoid's are
// right, so `J^R_X = R_Xᵀ J^L R_m` and `J^R_ω = J^L R_m` (docs/NUMERICS.md §21, *Jacobians, left*).
// `So3Spline` and `RdSpline` take integer nanoseconds alone: a record is answered by them only
// where `(u, h⁻¹)` is `(t / Δ, 10⁹ / Δ)` exactly, so that time quantization is never scored as
// basalt's error. Every conversion has a hand-computed test in `convert_test.cpp`.
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <sophus/se3.hpp>
#include <sophus/so3.hpp>

namespace oracle {

// The carrier of a record's values and of an answer's.
using Vec = std::vector<double>;
using Mat3 = Eigen::Matrix3d;

// A record basalt cannot take as it is: answered `nan` (docs/decisions/0072 decision 5).
struct Refused : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// A record whose `(u, h⁻¹)` is no integer-nanosecond time: a refusal of its own, counted apart.
struct OffGrid : Refused {
  using Refused::Refused;
};

// --- elements --------------------------------------------------------------------------------------

// Control `i` of `q` (helicoid `[w, x, y, z]` per control) as Sophus's `SO3` storage `[x, y, z, w]`.
inline std::array<double, 4> sophus_so3_params(const Vec& q, std::size_t i) {
  return {q.at(4 * i + 1), q.at(4 * i + 2), q.at(4 * i + 3), q.at(4 * i)};
}

// Control `i` of `(q, x)` as Sophus's `SE3` storage: the quaternion's, then the translation.
inline std::array<double, 7> sophus_se3_params(const Vec& q, const Vec& x, std::size_t i) {
  const auto r = sophus_so3_params(q, i);
  return {r[0], r[1], r[2], r[3], x.at(3 * i), x.at(3 * i + 1), x.at(3 * i + 2)};
}

inline Sophus::SO3d sophus_so3(const std::array<double, 4>& params) {
  return Sophus::SO3d(Eigen::Map<const Sophus::SO3d>(params.data()));
}

inline Vec sophus_to_helicoid_quat(const Eigen::Quaterniond& q) { return {q.w(), q.x(), q.y(), q.z()}; }

inline Vec sophus_to_helicoid_vec3(const Eigen::Vector3d& v) { return {v.x(), v.y(), v.z()}; }

// SE(3)'s tangent: helicoid `[ω; v]`, Sophus `[υ; ω]`.
inline Vec sophus_to_helicoid_twist(const Sophus::SE3d::Tangent& a) { return {a(3), a(4), a(5), a(0), a(1), a(2)}; }

// A 3 × 3 matrix as the corpus's array, column-major.
inline Vec column_major(const Mat3& m) { return Vec(m.data(), m.data() + 9); }

// --- Jacobian sides ------------------------------------------------------------------------------

// The value's right Jacobian in control `m` from basalt's left one: `Ad_{X⁻¹} J^L Ad_{X_m}`.
inline Mat3 right_value_jacobian(const Mat3& left, const Sophus::SO3d& value, const Sophus::SO3d& control) {
  return value.matrix().transpose() * left * control.matrix();
}

// A rate's (`ω`, `ω̇`: tangents with no output side) right Jacobian from basalt's left one:
// `J^L Ad_{X_m}`.
inline Mat3 right_rate_jacobian(const Mat3& left, const Sophus::SO3d& control) { return left * control.matrix(); }

// --- time ----------------------------------------------------------------------------------------

// What basalt's integer-nanosecond classes read: knot spacing `dt` and time `t` since the first
// knot, all `So3Spline`'s and `RdSpline`'s evaluations taking `u = double(t) / double(dt)` and
// `h⁻¹ = 1e9 / dt`.
struct Nanoseconds {
  std::int64_t dt, t;
};

// The `(dt, t)` at which basalt reads exactly `inv_h` and `u`, with `t / dt = u` exact in the
// reals; `OffGrid` where there is none.
inline Nanoseconds nanoseconds(double u, double inv_h) {
  const double dt = 1e9 / inv_h;
  const bool spacing = std::isfinite(dt) && dt >= 1 && dt < 0x1p53 && dt == std::nearbyint(dt) && 1e9 / dt == inv_h;
  const double t = u * dt;
  const bool time = spacing && t >= 0 && t < dt && t == std::nearbyint(t) && std::fma(u, dt, -t) == 0;
  if (!time) throw OffGrid("(u, h⁻¹) is no integer-nanosecond time");
  return {static_cast<std::int64_t>(dt), static_cast<std::int64_t>(t)};
}

// --- outputs -------------------------------------------------------------------------------------

// `values` scaled to unit size as the ids state: divided by `h⁻¹` once per derivative order, the
// shipped subject's `scaled`.
inline Vec scaled(Vec values, int order, double inv_h) {
  if (order == 0) return values;
  const double by = order == 1 ? inv_h : inv_h * inv_h;
  for (double& v : values) v /= by;
  return values;
}

}  // namespace oracle
