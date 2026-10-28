// The corpus's arrays to Sophus's types and back (docs/decisions/0072), at the scalar `T` a record
// is answered in: `double`, or `float` for an `@f32` record (docs/decisions/0077); and `eig3`'s and
// `svd3`'s to Eigen's, the image's own Eigen (docs/decisions/0083, 0101); `nearest_rotation`'s
// through Sophus's `SO3::fitToSO3` (docs/decisions/0102).
//
// helicoid is rotation-first and `w`-first (docs/decisions/0002); Sophus is translation-first:
// SE(2) `[υ; θ]`, SE(3) `[υ; ω]`, Sim(3) `[υ; ω; σ]`. Its Sim(3) stores `√s · q` (`RxSO3`), so the
// scale is read back as `log(‖√s · q‖²)`. Eigen's quaternion constructor takes `(w, x, y, z)` and
// stores `[x, y, z, w]`. The corpus's values are binary64: an input crosses into `T` by `narrow`,
// which refuses to round, and an answer back by the exact widening. Every conversion is a named
// function with a hand-computed test in `convert_test.cpp`, run when the image is built.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/SVD>
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

// --- eig3: Eigen's symmetric solver (docs/decisions/0083) -----------------------------------------

// A corpus `A` (N × N, column-major) as Eigen's `Matrix<T, N, N>`, column-major too;
// `SelfAdjointEigenSolver` reads its lower triangle, `JacobiSVD` all of it. N = 3 for `eig3` and
// `svd3`, 4, 9, 16 for `eig_sym_n*` (docs/decisions/0111).
template <class T, int N>
Eigen::Matrix<T, N, N> helicoid_to_eigen_matrix(const Vec& a) {
  if (a.size() != static_cast<std::size_t>(N * N)) {
    throw std::invalid_argument("an N × N matrix has N² values");
  }
  Eigen::Matrix<T, N, N> m;
  for (int k = 0; k < N * N; ++k) m.data()[k] = narrow<T>(a[static_cast<std::size_t>(k)]);
  return m;
}

template <class T>
Eigen::Matrix<T, 3, 3> helicoid_to_eigen_matrix3(const Vec& a) {
  return helicoid_to_eigen_matrix<T, 3>(a);
}

struct HelicoidEig3 {
  Vec lambda, v;
};

// A finished `SelfAdjointEigenSolver` as `eig3`'s (or `eig_sym_n*`'s) outputs: `lambda` ascending
// (Eigen sorts them; the stable sort states the contract, as nalgebra's runner does), `V` the
// matching columns, column-major, each column's sign Eigen's. No convergence (`info()`) is a
// refusal, so `nan`.
template <class T, int N>
HelicoidEig3 eigen_to_helicoid_eig(const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, N, N>>& es) {
  if (es.info() != Eigen::Success) throw Refused("SelfAdjointEigenSolver: no convergence");
  const auto& values = es.eigenvalues();
  const auto& vectors = es.eigenvectors();
  std::array<int, N> order{};
  for (int i = 0; i < N; ++i) order[static_cast<std::size_t>(i)] = i;
  std::stable_sort(order.begin(), order.end(), [&](int i, int j) { return values(i) < values(j); });
  HelicoidEig3 out;
  for (const int i : order) {
    out.lambda.push_back(static_cast<double>(values(i)));
    for (int r = 0; r < N; ++r) out.v.push_back(static_cast<double>(vectors(r, i)));
  }
  return out;
}

template <class T>
HelicoidEig3 eigen_to_helicoid_eig3(const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, 3, 3>>& es) {
  return eigen_to_helicoid_eig<T, 3>(es);
}

// --- svd3: Eigen's one-sided Jacobi SVD (docs/decisions/0101) ------------------------------------

struct HelicoidSvd3 {
  Vec u, sigma, v;
};

// A finished `JacobiSVD` (`ComputeFullU | ComputeFullV`) in record 0101 decision 1's convention:
// the triples sorted descending by `σ` (Eigen sorts; the stable sort states it); where `det V < 0`,
// `v₃` and `σ₃` negated; then where `det U < 0`, `u₃` and `σ₃`. So `U, V ∈ SO(3)` and
// `σ₁ ≥ σ₂ ≥ |σ₃|`. A non-finite input (`info()`) is a refusal, so `nan`.
template <class T>
HelicoidSvd3 eigen_to_helicoid_svd3(const Eigen::JacobiSVD<Eigen::Matrix<T, 3, 3>>& svd) {
  if (svd.info() != Eigen::Success) throw Refused("JacobiSVD: invalid input");
  const auto& values = svd.singularValues();
  std::array<int, 3> order{0, 1, 2};
  std::stable_sort(order.begin(), order.end(), [&](int i, int j) { return values(i) > values(j); });
  Eigen::Matrix<T, 3, 3> u, v;
  Eigen::Matrix<T, 3, 1> sigma;
  for (int k = 0; k < 3; ++k) {
    u.col(k) = svd.matrixU().col(order[k]);
    v.col(k) = svd.matrixV().col(order[k]);
    sigma(k) = values(order[k]);
  }
  for (Eigen::Matrix<T, 3, 3>* m : {&v, &u}) {
    if (m->determinant() < T(0)) {
      m->col(2) = -m->col(2);
      sigma(2) = -sigma(2);
    }
  }
  return {Vec(u.data(), u.data() + 9), {static_cast<double>(sigma(0)), static_cast<double>(sigma(1)),
                                        static_cast<double>(sigma(2))},
          Vec(v.data(), v.data() + 9)};
}

// --- nearest_rotation: Sophus's `fitToSO3` (docs/decisions/0102) --------------------------

// `SO3<T>::fitToSO3(M).matrix()`, column-major: Sophus's `makeRotationMatrix` (`JacobiSVD`, then
// `U diag(1, 1, det UVᵀ) Vᵀ`), stored as a quaternion by `SO3(R)`, whose orthogonality check is a
// refusal, and read back as a matrix.
template <class T>
Vec helicoid_nearest_rotation(const Vec& m) {
  return sophus_to_helicoid_matrix3<T>(Sophus::SO3<T>::fitToSO3(helicoid_to_eigen_matrix3<T>(m)).matrix());
}

}  // namespace oracle
