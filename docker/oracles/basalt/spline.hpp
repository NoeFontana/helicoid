// basalt-headers' answers to the spline ids (docs/decisions/0116 decision 4), one function per id,
// shared by `runner.cpp` and `convert_test.cpp`.
//
// Values and rates are `CeresSplineHelper<K>`'s, which takes `u` and `h⁻¹` as given and is generic
// over the Sophus group: `evaluate_lie` on `SO3` and on `SE3` (cumulative on SE(3); basalt's
// `Se3Spline` class is split and answers nothing here), `evaluate` on R³. Jacobians are the
// integer-nanosecond classes' (`So3Spline<K>`, `RdSpline<3, K>`), as is the Gram matrix; a record
// off their time grid is `OffGrid` (`convert.hpp`, `nanoseconds`).
#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <basalt/spline/ceres_spline_helper.h>
#include <basalt/spline/rd_spline.h>
#include <basalt/spline/so3_spline.h>

#include "convert.hpp"

namespace oracle {

using Fields = std::map<std::string, Vec>;

inline const Vec& take(const Fields& in, const char* key) {
  const auto it = in.find(key);
  if (it == in.end()) throw std::invalid_argument(std::string("no input `") + key + "`");
  return it->second;
}

inline double scalar(const Fields& in, const char* key) {
  const Vec& v = take(in, key);
  if (v.size() != 1) throw std::invalid_argument(std::string("`") + key + "` is not a scalar");
  return v[0];
}

// The window's order: `q`'s quaternions, or `p`'s points.
inline std::size_t order_of(const Fields& in, const char* key, std::size_t per) {
  const Vec& v = take(in, key);
  if (v.empty() || v.size() % per != 0) throw std::invalid_argument(std::string("`") + key + "` is no whole window");
  return v.size() / per;
}

template <int K>
using Order = std::integral_constant<int, K>;

// `f(Order<K>{})` at the record's `K`: the group orders of docs/decisions/0113.
template <class F>
Fields by_group_order(std::size_t k, F f) {
  switch (k) {
    case 4:
      return f(Order<4>{});
    case 5:
      return f(Order<5>{});
    case 6:
      return f(Order<6>{});
    default:
      throw std::invalid_argument("no group spline of order " + std::to_string(k));
  }
}

// The same over Rⁿ's orders.
template <class F>
Fields by_rn_order(std::size_t k, F f) {
  switch (k) {
    case 2:
      return f(Order<2>{});
    case 3:
      return f(Order<3>{});
    default:
      return by_group_order(k, f);
  }
}

// --- CeresSplineHelper: values and rates at `u` ---------------------------------------------------

// `evaluate_lie` on SO(3) at derivative order `D`.
template <int D>
Fields so3_rate(const Fields& in) {
  const Vec& q = take(in, "q");
  const double u = scalar(in, "u"), inv_h = scalar(in, "inv_h");
  return by_group_order(order_of(in, "q", 4), [&](auto n) {
    constexpr int K = decltype(n)::value;
    std::array<std::array<double, 4>, K> knots;
    std::array<const double*, K> at{};
    for (int i = 0; i < K; ++i) {
      knots[i] = sophus_so3_params(q, i);
      at[i] = knots[i].data();
    }
    Sophus::SO3d x;
    Eigen::Vector3d rate;
    using Helper = basalt::CeresSplineHelper<K>;
    if constexpr (D == 0) {
      Helper::template evaluate_lie<double, Sophus::SO3>(at.data(), u, inv_h, &x);
      return Fields{{"q", sophus_to_helicoid_quat(x.unit_quaternion())}};
    } else if constexpr (D == 1) {
      Helper::template evaluate_lie<double, Sophus::SO3>(at.data(), u, inv_h, nullptr, &rate);
      return Fields{{"omega", scaled(sophus_to_helicoid_vec3(rate), 1, inv_h)}};
    } else {
      Helper::template evaluate_lie<double, Sophus::SO3>(at.data(), u, inv_h, nullptr, nullptr, &rate);
      return Fields{{"alpha", scaled(sophus_to_helicoid_vec3(rate), 2, inv_h)}};
    }
  });
}

// `evaluate_lie` on SE(3): the cumulative spline on the group, its twists `[υ; ω]` read back
// rotation-first.
template <int D>
Fields se3_rate(const Fields& in) {
  const Vec &q = take(in, "q"), &t = take(in, "x");
  const double u = scalar(in, "u"), inv_h = scalar(in, "inv_h");
  const std::size_t k = order_of(in, "q", 4);
  if (t.size() != 3 * k) throw std::invalid_argument("`x` is not one translation per control");
  return by_group_order(k, [&](auto n) {
    constexpr int K = decltype(n)::value;
    std::array<std::array<double, 7>, K> knots;
    std::array<const double*, K> at{};
    for (int i = 0; i < K; ++i) {
      knots[i] = sophus_se3_params(q, t, i);
      at[i] = knots[i].data();
    }
    Sophus::SE3d x;
    Sophus::SE3d::Tangent rate;
    using Helper = basalt::CeresSplineHelper<K>;
    if constexpr (D == 0) {
      Helper::template evaluate_lie<double, Sophus::SE3>(at.data(), u, inv_h, &x);
      return Fields{{"q", sophus_to_helicoid_quat(x.unit_quaternion())},
                    {"x", sophus_to_helicoid_vec3(x.translation())}};
    } else if constexpr (D == 1) {
      Helper::template evaluate_lie<double, Sophus::SE3>(at.data(), u, inv_h, nullptr, &rate);
      return Fields{{"twist", scaled(sophus_to_helicoid_twist(rate), 1, inv_h)}};
    } else {
      Helper::template evaluate_lie<double, Sophus::SE3>(at.data(), u, inv_h, nullptr, nullptr, &rate);
      return Fields{{"twist", scaled(sophus_to_helicoid_twist(rate), 2, inv_h)}};
    }
  });
}

// `evaluate<double, 3, D>`: the non-cumulative sum over the window.
template <int D>
Fields rn_rate(const Fields& in) {
  const Vec& p = take(in, "p");
  const double u = scalar(in, "u"), inv_h = scalar(in, "inv_h");
  return by_rn_order(order_of(in, "p", 3), [&](auto n) {
    constexpr int K = decltype(n)::value;
    std::array<const double*, K> at{};
    for (int i = 0; i < K; ++i) at[i] = p.data() + 3 * i;
    Eigen::Vector3d out;
    basalt::CeresSplineHelper<K>::template evaluate<double, 3, D>(at.data(), u, inv_h, &out);
    const char* field = D == 0 ? "p" : D == 1 ? "v" : "a";
    return Fields{{field, scaled(sophus_to_helicoid_vec3(out), D, inv_h)}};
  });
}

// --- So3Spline, RdSpline: Jacobians on the integer-nanosecond grid --------------------------------

// `So3Spline<K>`'s Jacobian of the value (`D = 0`), `ω` or `ω̇` in every control, made right-sided.
template <int D>
Fields so3_jacobian(const Fields& in) {
  const Vec& q = take(in, "q");
  const double inv_h = scalar(in, "inv_h");
  const Nanoseconds ns = nanoseconds(scalar(in, "u"), inv_h);
  return by_group_order(order_of(in, "q", 4), [&](auto n) {
    constexpr int K = decltype(n)::value;
    basalt::So3Spline<K> spline(ns.dt);
    for (int i = 0; i < K; ++i) spline.knotsPushBack(sophus_so3(sophus_so3_params(q, i)));
    typename basalt::So3Spline<K>::JacobianStruct j;
    Sophus::SO3d x;
    if constexpr (D == 0) {
      x = spline.evaluate(ns.t, &j);
    } else if constexpr (D == 1) {
      spline.velocityBody(ns.t, &j);
    } else {
      spline.accelerationBody(ns.t, &j);
    }
    if (j.start_idx != 0) throw std::logic_error("So3Spline evaluated another segment");
    Vec jac;
    for (int m = 0; m < K; ++m) {
      const Sophus::SO3d& control = spline.getKnot(m);
      const Mat3 right = D == 0 ? right_value_jacobian(j.d_val_d_knot[m], x, control)
                                : right_rate_jacobian(j.d_val_d_knot[m], control);
      const Vec entries = column_major(right);
      jac.insert(jac.end(), entries.begin(), entries.end());
    }
    return Fields{{"jac", scaled(jac, D, inv_h)}};
  });
}

// `RdSpline<3, K>`'s `d_val_d_knot`: the `k` of each control's `k I`.
template <int D>
Fields rn_jacobian(const Fields& in) {
  const Vec& p = take(in, "p");
  const double inv_h = scalar(in, "inv_h");
  const Nanoseconds ns = nanoseconds(scalar(in, "u"), inv_h);
  return by_rn_order(order_of(in, "p", 3), [&](auto n) {
    constexpr int K = decltype(n)::value;
    if constexpr (K == 2 && D == 2) {
      // `evaluate<2>` reads `pow_inv_dt_[2]` of a `std::array` of 2.
      throw Refused("RdSpline<3, 2> has no acceleration");
      return Fields{};
    } else {
      basalt::RdSpline<3, K> spline(ns.dt);
      for (int i = 0; i < K; ++i) spline.knotsPushBack(Eigen::Vector3d(p[3 * i], p[3 * i + 1], p[3 * i + 2]));
      typename basalt::RdSpline<3, K>::JacobianStruct j;
      spline.template evaluate<D>(ns.t, &j);
      if (j.start_idx != 0) throw std::logic_error("RdSpline evaluated another segment");
      return Fields{{"k", scaled(Vec(j.d_val_d_knot.begin(), j.d_val_d_knot.end()), D, inv_h)}};
    }
  });
}

// `G⁽ʳ⁾`, column-major: `RdSpline<1, K>::evaluateIntegralSquared<r>` is `½ cᵀ Q c` with its
// gradient `Q c`, so `G = Q / 2`'s column `t` is half the gradient at the unit knot `e_t`.
inline Fields rn_gram(const Fields& in) {
  const double k = scalar(in, "k"), r = scalar(in, "r");
  const Nanoseconds ns = nanoseconds(0, scalar(in, "inv_h"));
  return by_rn_order(static_cast<std::size_t>(k), [&](auto n) {
    constexpr int K = decltype(n)::value;
    const auto column = [&](int t, auto order) {
      constexpr int R = decltype(order)::value;
      basalt::RdSpline<1, K> spline(ns.dt);
      for (int i = 0; i < K; ++i) spline.knotsPushBack(Eigen::Matrix<double, 1, 1>(i == t ? 1.0 : 0.0));
      typename basalt::RdSpline<1, K>::IntegralJacobianStruct j;
      spline.template evaluateIntegralSquared<R>(&j);
      Vec half;
      for (int s = 0; s < K; ++s) half.push_back(0.5 * j.d_int_d_knot.at(s)(0));
      return half;
    };
    Vec g;
    for (int t = 0; t < K; ++t) {
      Vec c;
      if (r == 1) {
        c = column(t, Order<1>{});
      } else if (r == 2) {
        c = column(t, Order<2>{});
      } else {
        throw std::invalid_argument("no Gram matrix of derivative order " + std::to_string(r));
      }
      g.insert(g.end(), c.begin(), c.end());
    }
    return Fields{{"g", g}};
  });
}

}  // namespace oracle
