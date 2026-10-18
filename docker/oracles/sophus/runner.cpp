// The Sophus C++ oracle (docs/PHASE1.md §7, docs/decisions/0072): reads corpus JSONL and writes one
// answer file per function id it answers. It never reads a reference and never scores;
// `cargo xtask conformance --oracle sophus` does, with the harness's exact metric.
//
//   sophus_runner --version                               the pin, one line: the result rows' subject_version
//   sophus_runner [--precision f32] --out DIR FILE.jsonl…  DIR/<id>.jsonl for each FILE whose id is answered
//
// An answer line is `{"id":N,"out":{"<field>":["<hex float>",…]}}`, fields flat in the reference's
// order (matrices column-major). In binary64 every record is answered with `Sophus::*<double>`;
// with `--precision f32` the `@f32` records alone, with `Sophus::*<float>` (docs/decisions/0077).
// `eig3` is the image's Eigen, `SelfAdjointEigenSolver<Matrix3<T>>` (docs/decisions/0083); `svd3`
// its `JacobiSVD<Matrix3<T>>` (docs/decisions/0101); `nearest_rotation` `SO3<T>::fitToSO3`
// (docs/decisions/0102, draft).
// A record Sophus refuses (a failed `SOPHUS_ENSURE`, or Eigen's solver not converging) is answered
// with `nan` in every field: answered, so no stratum is the runner's choice, and non-finite, so it
// never sets a bar.
#include <gnu/libc-version.h>

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>
#include <sophus/interpolate.hpp>

#include "convert.hpp"
#include "hexfloat.hpp"

#ifndef ORACLE_SOPHUS
#error "ORACLE_SOPHUS (the pinned Sophus version and commit) is set by the Dockerfile"
#endif

namespace oracle {
namespace {

using Fields = std::map<std::string, Vec>;

// --- the ids -------------------------------------------------------------------------------------

const Vec& take(const Fields& in, const std::string& key, std::size_t len) {
  const auto it = in.find(key);
  if (it == in.end()) throw std::invalid_argument("no input `" + key + "`");
  if (it->second.size() != len) {
    throw std::invalid_argument("`" + key + "` has " + std::to_string(it->second.size()) +
                                " values, expected " + std::to_string(len));
  }
  return it->second;
}

double scalar(const Fields& in, const std::string& key) { return take(in, key, 1)[0]; }

// Which Jacobian an id asks for. Sophus has the left one and its inverse; the right ones are
// answered at the negated tangent, `J_r(τ) = J_l(−τ)` (docs/NUMERICS.md §1).
enum class Jacobian { Left, Right, LeftInverse, RightInverse };

template <class V>
V jacobian_tangent(Jacobian which, const V& tau) {
  return which == Jacobian::Right || which == Jacobian::RightInverse ? V(-tau) : tau;
}

bool is_inverse(Jacobian which) {
  return which == Jacobian::LeftInverse || which == Jacobian::RightInverse;
}

// The ids at the scalar `T`: every input narrowed to `T` (`convert.hpp`), every operation Sophus's
// at `T`.
template <class T>
struct At {
  using SO2 = Sophus::SO2<T>;
  using SE2 = Sophus::SE2<T>;
  using SO3 = Sophus::SO3<T>;
  using SE3 = Sophus::SE3<T>;
  using Sim3 = Sophus::Sim3<T>;

  static T number(const Fields& in, const std::string& key) { return narrow<T>(scalar(in, key)); }

  static Fields so3_jacobian(const Fields& in, Jacobian which) {
    const Eigen::Matrix<T, 3, 1> phi = jacobian_tangent(which, helicoid_to_sophus_vec3<T>(take(in, "phi", 3)));
    const Eigen::Matrix<T, 3, 3> j = is_inverse(which) ? SO3::leftJacobianInverse(phi) : SO3::leftJacobian(phi);
    return {{"J", sophus_to_helicoid_matrix3<T>(j)}};
  }

  static Fields se3_jacobian(const Fields& in, Jacobian which) {
    const typename SE3::Tangent tau =
        jacobian_tangent(which, helicoid_to_sophus_tangent<T>(take(in, "tau", 6), kSe3));
    const typename SE3::Adjoint j = is_inverse(which) ? SE3::leftJacobianInverse(tau) : SE3::leftJacobian(tau);
    return {{"J", sophus_to_helicoid_matrix<T>(j, kSe3)}};
  }

  static Fields se2_out(const SE2& g) {
    const Eigen::Matrix<T, 2, 1> t = g.translation();
    return {{"t", {static_cast<double>(t.x()), static_cast<double>(t.y())}}, {"z", sophus_to_helicoid_so2<T>(g.so2())}};
  }

  static Fields se3_out(const SE3& g) {
    return {{"q", sophus_to_helicoid_quat<T>(g.unit_quaternion())}, {"x", sophus_to_helicoid_vec3<T>(g.translation())}};
  }

  static SE2 se2_in(const Fields& in) { return helicoid_to_sophus_se2<T>(take(in, "z", 2), take(in, "t", 2)); }

  static SE3 se3_in(const Fields& in) { return helicoid_to_sophus_se3<T>(take(in, "q", 4), take(in, "x", 3)); }

  static Sim3 sim3_in(const Fields& in) {
    return helicoid_to_sophus_sim3<T>(take(in, "q", 4), take(in, "x", 3), scalar(in, "sigma"));
  }

  static Fields evaluate(const std::string& id, const Fields& in) {
    using J = Jacobian;
    if (id == "so2_exp") return {{"z", sophus_to_helicoid_so2<T>(SO2::exp(number(in, "theta")))}};
    if (id == "so2_log") return {{"theta", {static_cast<double>(helicoid_to_sophus_so2<T>(take(in, "z", 2)).log())}}};
    if (id == "se2_exp") return se2_out(SE2::exp(helicoid_to_sophus_tangent<T>(take(in, "tau", 3), kSe2)));
    if (id == "se2_log") return {{"tau", sophus_to_helicoid_tangent<T>(se2_in(in).log(), kSe2)}};
    if (id == "se2_ad") return {{"Ad", sophus_to_helicoid_matrix<T>(se2_in(in).Adj(), kSe2)}};
    if (id == "so3_exp") {
      return {{"q", sophus_to_helicoid_quat<T>(SO3::exp(helicoid_to_sophus_vec3<T>(take(in, "phi", 3))).unit_quaternion())}};
    }
    if (id == "so3_log") return {{"phi", sophus_to_helicoid_vec3<T>(helicoid_to_sophus_so3<T>(take(in, "q", 4)).log())}};
    if (id == "so3_act") {
      return {{"Rp", sophus_to_helicoid_vec3<T>(helicoid_to_sophus_so3<T>(take(in, "q", 4)) *
                                                 helicoid_to_sophus_vec3<T>(take(in, "p", 3)))}};
    }
    // `SO3(R)`: Eigen's quaternion of the matrix, after Sophus's orthogonality check.
    if (id == "so3_from_matrix") {
      const Vec r = take(in, "R", 9);
      Eigen::Matrix<T, 3, 3> m;
      for (int i = 0; i < 9; ++i) m.data()[i] = narrow<T>(r[static_cast<std::size_t>(i)]);
      return {{"q", sophus_to_helicoid_quat<T>(SO3(m).unit_quaternion())}};
    }
    if (id == "so3_jl") return so3_jacobian(in, J::Left);
    if (id == "so3_jr") return so3_jacobian(in, J::Right);
    if (id == "so3_jl_inv") return so3_jacobian(in, J::LeftInverse);
    if (id == "so3_jr_inv") return so3_jacobian(in, J::RightInverse);
    // `γ(X₀, X₁, t) = X₀ Exp(t Log(X₀⁻¹ X₁))` (docs/NUMERICS.md §10) is `Sophus::interpolate`.
    if (id == "so3_geodesic") {
      const SO3 q = Sophus::interpolate(helicoid_to_sophus_so3<T>(take(in, "q0", 4)),
                                        helicoid_to_sophus_so3<T>(take(in, "q1", 4)), number(in, "t"));
      return {{"q", sophus_to_helicoid_quat<T>(q.unit_quaternion())}};
    }
    if (id == "se3_geodesic") {
      return se3_out(Sophus::interpolate(helicoid_to_sophus_se3<T>(take(in, "q0", 4), take(in, "x0", 3)),
                                         helicoid_to_sophus_se3<T>(take(in, "q1", 4), take(in, "x1", 3)),
                                         number(in, "t")));
    }
    if (id == "sen3_exp_n1") return se3_out(SE3::exp(helicoid_to_sophus_tangent<T>(take(in, "tau", 6), kSe3)));
    if (id == "sen3_log_n1") return {{"tau", sophus_to_helicoid_tangent<T>(se3_in(in).log(), kSe3)}};
    if (id == "sen3_ad_n1") return {{"Ad", sophus_to_helicoid_matrix<T>(se3_in(in).Adj(), kSe3)}};
    if (id == "sen3_jl_n1") return se3_jacobian(in, J::Left);
    if (id == "sen3_jr_n1") return se3_jacobian(in, J::Right);
    if (id == "sen3_jl_inv_n1") return se3_jacobian(in, J::LeftInverse);
    if (id == "sen3_jr_inv_n1") return se3_jacobian(in, J::RightInverse);
    if (id == "sim3_exp") {
      const auto x = sophus_to_helicoid_sim3<T>(Sim3::exp(helicoid_to_sophus_tangent<T>(take(in, "tau", 7), kSim3)));
      return {{"q", x.q}, {"x", x.x}, {"sigma", {x.sigma}}};
    }
    if (id == "sim3_log") return {{"tau", sophus_to_helicoid_tangent<T>(sim3_in(in).log(), kSim3)}};
    if (id == "sim3_ad") return {{"Ad", sophus_to_helicoid_matrix<T>(sim3_in(in).Adj(), kSim3)}};
    if (id == "sim3_act") {
      return {{"Xp", sophus_to_helicoid_vec3<T>(sim3_in(in) * helicoid_to_sophus_vec3<T>(take(in, "p", 3)))}};
    }
    // Not Sophus's: the image's Eigen, the solver a Sophus user has (docs/decisions/0083).
    if (id == "eig3") {
      const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, 3, 3>> es(helicoid_to_eigen_matrix3<T>(take(in, "A", 9)));
      const HelicoidEig3 e = eigen_to_helicoid_eig3<T>(es);
      return {{"V", e.v}, {"lambda", e.lambda}};
    }
    // The image's Eigen too, signed by the runner (docs/decisions/0101).
    if (id == "svd3") {
      const Eigen::JacobiSVD<Eigen::Matrix<T, 3, 3>> svd(helicoid_to_eigen_matrix3<T>(take(in, "A", 9)),
                                                        Eigen::ComputeFullU | Eigen::ComputeFullV);
      const HelicoidSvd3 s = eigen_to_helicoid_svd3<T>(svd);
      return {{"U", s.u}, {"sigma", s.sigma}, {"V", s.v}};
    }
    if (id == "nearest_rotation") return {{"R", helicoid_nearest_rotation<T>(take(in, "M", 9))}};
    throw std::logic_error("no evaluation of `" + id + "`");
  }
};

// One id: the fields and lengths it answers, which a refusal fills with `nan`, and whether it is
// answered in binary32: the ids whose corpus file has `@f32` strata (docs/decisions/0077).
struct Id {
  std::string name;
  std::vector<std::pair<std::string, std::size_t>> out;
  bool f32;
};

const std::vector<Id>& table() {
  static const std::vector<Id> ids = {
      {"so2_exp", {{"z", 2}}, false},
      {"so2_log", {{"theta", 1}}, false},
      {"se2_exp", {{"t", 2}, {"z", 2}}, false},
      {"se2_log", {{"tau", 3}}, false},
      {"se2_ad", {{"Ad", 9}}, false},
      {"so3_exp", {{"q", 4}}, false},
      {"so3_log", {{"phi", 3}}, false},
      {"so3_act", {{"Rp", 3}}, false},
      {"so3_from_matrix", {{"q", 4}}, true},
      {"so3_jl", {{"J", 9}}, false},
      {"so3_jr", {{"J", 9}}, false},
      {"so3_jl_inv", {{"J", 9}}, false},
      {"so3_jr_inv", {{"J", 9}}, false},
      {"so3_geodesic", {{"q", 4}}, false},
      {"se3_geodesic", {{"q", 4}, {"x", 3}}, false},
      {"sen3_exp_n1", {{"q", 4}, {"x", 3}}, false},
      {"sen3_log_n1", {{"tau", 6}}, false},
      {"sen3_ad_n1", {{"Ad", 36}}, false},
      {"sen3_jl_n1", {{"J", 36}}, false},
      {"sen3_jr_n1", {{"J", 36}}, false},
      {"sen3_jl_inv_n1", {{"J", 36}}, false},
      {"sen3_jr_inv_n1", {{"J", 36}}, false},
      {"sim3_exp", {{"q", 4}, {"x", 3}, {"sigma", 1}}, true},
      {"sim3_log", {{"tau", 7}}, true},
      {"sim3_ad", {{"Ad", 49}}, true},
      {"sim3_act", {{"Xp", 3}}, true},
      {"eig3", {{"V", 9}, {"lambda", 3}}, true},
      {"svd3", {{"U", 9}, {"sigma", 3}, {"V", 9}}, true},
      {"nearest_rotation", {{"R", 9}}, true},
  };
  return ids;
}

// The id `name` answered at `f32` or not, or none.
const Id* find(const std::string& name, bool f32) {
  for (const Id& id : table()) {
    if (id.name == name && (id.f32 || !f32)) return &id;
  }
  return nullptr;
}

// --- the files -----------------------------------------------------------------------------------

// A record's `in` object: each entry a hex-float string or an array of them; a matrix's `shape`
// is the corpus's and not an input.
Fields inputs(const nlohmann::json& object) {
  Fields fields;
  for (const auto& [key, value] : object.items()) {
    if (key == "shape") continue;
    Vec values;
    if (value.is_string()) {
      values.push_back(parse_hex(value.get<std::string>()));
    } else if (value.is_array()) {
      for (const auto& v : value) {
        if (!v.is_string()) throw std::invalid_argument("`" + key + "`: not a string");
        values.push_back(parse_hex(v.get<std::string>()));
      }
    } else {
      throw std::invalid_argument("`" + key + "`: neither a string nor an array");
    }
    fields.emplace(key, std::move(values));
  }
  return fields;
}

// The answer to one record at `float` or `double`: Sophus's, its shape checked, or `nan`
// throughout if Sophus refuses.
Fields answer(const Id& id, const Fields& in, bool f32, std::size_t& refused) {
  Fields out;
  try {
    out = f32 ? At<float>::evaluate(id.name, in) : At<double>::evaluate(id.name, in);
  } catch (const Refused&) {
    ++refused;
    for (const auto& [field, len] : id.out) out[field] = Vec(len, std::nan(""));
  }
  if (out.size() != id.out.size()) throw std::logic_error(id.name + ": wrong output fields");
  for (const auto& [field, len] : id.out) {
    const auto it = out.find(field);
    if (it == out.end() || it->second.size() != len) {
      throw std::logic_error(id.name + ": output `" + field + "` has the wrong length");
    }
  }
  return out;
}

std::string line_of(std::uint64_t record, const Fields& out) {
  std::ostringstream line;
  line << "{\"id\":" << record << ",\"out\":{";
  bool first_field = true;
  for (const auto& [field, values] : out) {
    line << (first_field ? "" : ",") << '"' << field << "\":[";
    for (std::size_t i = 0; i < values.size(); ++i) line << (i ? "," : "") << '"' << format_hex(values[i]) << '"';
    line << ']';
    first_field = false;
  }
  line << "}}\n";
  return line.str();
}

bool is_f32_stratum(const std::string& stratum) {
  const std::string suffix = "@f32";
  return stratum.size() >= suffix.size() && stratum.compare(stratum.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// Answers every line of `path` at its precision (in binary32 the `@f32` lines alone) into
// `target`; returns the number of refused records.
std::size_t answer_file(const Id& id, const std::string& path, const std::string& target, bool f32) {
  std::ifstream file(path);
  if (!file) throw std::runtime_error(path + ": cannot read");
  std::ofstream out(target);
  if (!out) throw std::runtime_error(target + ": cannot write");
  std::string text;
  std::size_t number = 0, refused = 0;
  while (std::getline(file, text)) {
    ++number;
    try {
      const auto record = nlohmann::json::parse(text);
      if (f32 && !is_f32_stratum(record.at("stratum").get<std::string>())) continue;
      const std::uint64_t record_id = record.at("id").get<std::uint64_t>();
      out << line_of(record_id, answer(id, inputs(record.at("in")), f32, refused));
    } catch (const std::exception& e) {
      throw std::runtime_error(path + ": line " + std::to_string(number) + ": " + e.what());
    }
  }
  if (!out.flush()) throw std::runtime_error(target + ": write failed");
  return refused;
}

std::string pin() {
  std::ostringstream text;
  text << "sophus@" << ORACLE_SOPHUS << " eigen@" << EIGEN_WORLD_VERSION << '.' << EIGEN_MAJOR_VERSION << '.'
       << EIGEN_MINOR_VERSION << " gcc@" << __GNUC__ << '.' << __GNUC_MINOR__ << '.' << __GNUC_PATCHLEVEL__
       << " glibc@" << gnu_get_libc_version();
  return text.str();
}

int run(std::vector<std::string> args) {
  if (args.size() == 1 && args[0] == "--version") {
    std::cout << pin() << '\n';
    return std::cout.flush() ? 0 : 1;
  }
  const bool f32 = args.size() >= 2 && args[0] == "--precision" && args[1] == "f32";
  if (f32) args.erase(args.begin(), args.begin() + 2);
  if (args.size() < 3 || args[0] != "--out") {
    std::cerr << "usage: sophus_runner --version | [--precision f32] --out DIR FILE.jsonl…\n";
    return 2;
  }
  const std::string dir = args[1];
  std::vector<std::string> answered;
  std::size_t skipped = 0, refused = 0;
  for (std::size_t i = 2; i < args.size(); ++i) {
    const std::string& path = args[i];
    const std::size_t slash = path.find_last_of('/');
    std::string stem = path.substr(slash == std::string::npos ? 0 : slash + 1);
    const std::size_t dot = stem.rfind(".jsonl");
    if (dot == std::string::npos || dot + 6 != stem.size()) throw std::invalid_argument(path + ": not a .jsonl file");
    stem.resize(dot);
    const Id* id = find(stem, f32);
    if (id == nullptr) {
      ++skipped;
      continue;
    }
    refused += answer_file(*id, path, dir + "/" + stem + ".jsonl", f32);
    answered.push_back(stem);
  }
  std::cerr << "answered " << answered.size() << " ids in " << (f32 ? "binary32" : "binary64") << "; " << refused
            << " records refused (nan); " << skipped << " files of other ids skipped\n";
  return 0;
}

}  // namespace
}  // namespace oracle

int main(int argc, char** argv) {
  try {
    return oracle::run(std::vector<std::string>(argv + 1, argv + argc));
  } catch (const std::exception& e) {
    std::cerr << "sophus_runner: " << e.what() << '\n';
    return 1;
  }
}
