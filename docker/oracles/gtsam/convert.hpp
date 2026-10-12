// The corpus's arrays to GTSAM's types and back (docs/decisions/0074).
//
// helicoid is rotation-first and `w`-first (docs/decisions/0002). GTSAM is rotation-first too for
// SO(3), SE_N(3) (`[ω; ρ₁; …; ρ_N]`, `Pose3` and `ExtendedPose3<K>`) and Sim(3) (`[ω; ρ; λ]`), so
// only SE(2) (`[v; ω]`) is permuted. Its `Similarity3` stores `(R, t, s)` acting as `s(Rp + t)`, so
// helicoid's translation is `s · t`. Its `Rot3` stores the matrix of the quaternion it is built
// from, unnormalized. Every conversion is a named function with a hand-computed test in
// `convert_test.cpp`, run when the image is built.
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtsam/geometry/ExtendedPose3.h>
#include <gtsam/geometry/Pose2.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/Rot2.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/geometry/Similarity3.h>

namespace oracle {

using Vec = std::vector<double>;

// --- quaternions and rotations -------------------------------------------------------------------

// `Rot3::Quaternion(w, x, y, z)` is the matrix `Eigen::Quaterniond::toRotationMatrix` forms, which
// assumes a unit quaternion and does not normalize.
inline gtsam::Rot3 helicoid_to_gtsam_rot3(const Vec& q) {
  return gtsam::Rot3::Quaternion(q.at(0), q.at(1), q.at(2), q.at(3));
}

// `toQuaternion` is Eigen's quaternion of the stored matrix; Eigen stores `[x, y, z, w]`.
inline Vec gtsam_to_helicoid_quat(const gtsam::Rot3& r) {
  const gtsam::Quaternion q = r.toQuaternion();
  return {q.w(), q.x(), q.y(), q.z()};
}

inline Eigen::Vector3d helicoid_to_gtsam_vec3(const Vec& v) {
  return Eigen::Vector3d(v.at(0), v.at(1), v.at(2));
}

inline Vec gtsam_to_helicoid_vec3(const Eigen::Vector3d& v) { return {v.x(), v.y(), v.z()}; }

// A matrix in helicoid's order as the corpus's array: column-major, as Eigen stores it.
template <int R, int C>
Vec gtsam_to_helicoid_matrix(const Eigen::Matrix<double, R, C>& m) {
  return Vec(m.data(), m.data() + m.size());
}

// The corpus's column-major `R` as GTSAM's matrix.
inline Eigen::Matrix3d helicoid_to_gtsam_matrix3(const Vec& r) {
  if (r.size() != 9) throw std::invalid_argument("a 3 × 3 matrix has 9 values");
  return Eigen::Map<const Eigen::Matrix3d>(r.data());
}

// --- tangents ------------------------------------------------------------------------------------

// A tangent in the same order in both libraries.
template <int N>
Eigen::Matrix<double, N, 1> helicoid_to_gtsam_tangent(const Vec& tau) {
  if (tau.size() != static_cast<std::size_t>(N)) {
    throw std::invalid_argument("a tangent has " + std::to_string(N) + " values");
  }
  return Eigen::Map<const Eigen::Matrix<double, N, 1>>(tau.data());
}

template <int N>
Vec gtsam_to_helicoid_tangent(const Eigen::Matrix<double, N, 1>& a) {
  return Vec(a.data(), a.data() + N);
}

// SE(2): helicoid `[θ; ρ₀; ρ₁]` is GTSAM `[v₀; v₁; ω]`; helicoid's component `i` sits at
// `kSe2At[i]`, and a Jacobian or adjoint's helicoid `(i, j)` is GTSAM's `(kSe2At[i], kSe2At[j])`.
inline constexpr std::array<int, 3> kSe2At{2, 0, 1};

inline Eigen::Vector3d helicoid_to_gtsam_se2_tangent(const Vec& tau) {
  if (tau.size() != 3) throw std::invalid_argument("an SE(2) tangent has 3 values");
  Eigen::Vector3d a;
  for (int i = 0; i < 3; ++i) a(kSe2At[i]) = tau[i];
  return a;
}

inline Vec gtsam_to_helicoid_se2_tangent(const Eigen::Vector3d& a) {
  Vec tau(3);
  for (int i = 0; i < 3; ++i) tau[i] = a(kSe2At[i]);
  return tau;
}

inline Vec gtsam_to_helicoid_se2_matrix(const Eigen::Matrix3d& m) {
  Vec out;
  out.reserve(9);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) out.push_back(m(kSe2At[i], kSe2At[j]));
  }
  return out;
}

// --- elements ------------------------------------------------------------------------------------

// helicoid's `z = [cos θ, sin θ]` through `Rot2::fromCosSin`, which normalizes a pair whose squared
// norm is off 1 by more than 10⁻¹⁰.
inline gtsam::Rot2 helicoid_to_gtsam_rot2(const Vec& z) { return gtsam::Rot2::fromCosSin(z.at(0), z.at(1)); }

inline Vec gtsam_to_helicoid_rot2(const gtsam::Rot2& r) { return {r.c(), r.s()}; }

inline gtsam::Pose2 helicoid_to_gtsam_pose2(const Vec& z, const Vec& t) {
  return gtsam::Pose2(helicoid_to_gtsam_rot2(z), gtsam::Point2(t.at(0), t.at(1)));
}

inline gtsam::Pose3 helicoid_to_gtsam_pose3(const Vec& q, const Vec& x) {
  return gtsam::Pose3(helicoid_to_gtsam_rot3(q), helicoid_to_gtsam_vec3(x));
}

// SE_N(3) for N = 2, 3: helicoid's `x` is `[x₁; …; x_N]`, the columns of GTSAM's `3 × K` block in
// order. (`NavState` is `ExtendedPose3<2>` with columns `(t, v)`; helicoid's `SE23` names its
// columns `(v, p)`. The group does not depend on the names, so the columns map in order.)
template <int K>
gtsam::ExtendedPose3<K> helicoid_to_gtsam_extended_pose3(const Vec& q, const Vec& x) {
  if (x.size() != static_cast<std::size_t>(3 * K)) {
    throw std::invalid_argument("an SE_N(3) translation block has " + std::to_string(3 * K) + " values");
  }
  const Eigen::Matrix<double, 3, K> block = Eigen::Map<const Eigen::Matrix<double, 3, K>>(x.data());
  return gtsam::ExtendedPose3<K>(helicoid_to_gtsam_rot3(q), block);
}

// `(q, x, σ)` as `Similarity3(R, x / e^σ, e^σ)`: GTSAM's translation is applied before the scale.
inline gtsam::Similarity3 helicoid_to_gtsam_similarity3(const Vec& q, const Vec& x, double sigma) {
  const double s = std::exp(sigma);
  return gtsam::Similarity3(helicoid_to_gtsam_rot3(q), helicoid_to_gtsam_vec3(x) / s, s);
}

struct HelicoidSim3 {
  Vec q, x;
  double sigma;
};

// Back: `x = s · t`, `σ = log s`.
inline HelicoidSim3 gtsam_to_helicoid_sim3(const gtsam::Similarity3& g) {
  const double s = g.scale();
  return {gtsam_to_helicoid_quat(g.rotation()), gtsam_to_helicoid_vec3(s * g.translation()), std::log(s)};
}

}  // namespace oracle
