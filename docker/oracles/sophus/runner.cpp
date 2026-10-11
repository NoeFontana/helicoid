// The Sophus C++ oracle (docs/PHASE1.md §7, docs/decisions/0072): reads corpus JSONL and writes one
// answer file per function id it answers. It never reads a reference and never scores;
// `cargo xtask conformance --oracle sophus` does, with the harness's exact metric.
//
//   sophus_runner --version               the pin, one line: the result rows' subject_version
//   sophus_runner --out DIR FILE.jsonl…   DIR/<id>.jsonl for each FILE whose id is answered
//
// An answer line is `{"id":N,"out":{"<field>":["<hex float>",…]}}`, fields flat in the reference's
// order (matrices column-major). A record Sophus refuses (a failed `SOPHUS_ENSURE`) is answered with
// `nan` in every field: answered, so no stratum is the runner's choice, and non-finite, so it never
// sets a bar.
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

template <class T>
T jacobian_tangent(Jacobian which, const T& tau) {
  return which == Jacobian::Right || which == Jacobian::RightInverse ? T(-tau) : tau;
}

bool is_inverse(Jacobian which) {
  return which == Jacobian::LeftInverse || which == Jacobian::RightInverse;
}

Fields so3_jacobian(const Fields& in, Jacobian which) {
  const Eigen::Vector3d phi = jacobian_tangent(which, helicoid_to_sophus_vec3(take(in, "phi", 3)));
  const Eigen::Matrix3d j = is_inverse(which) ? Sophus::SO3d::leftJacobianInverse(phi)
                                              : Sophus::SO3d::leftJacobian(phi);
  return {{"J", sophus_to_helicoid_matrix3(j)}};
}

Fields se3_jacobian(const Fields& in, Jacobian which) {
  const Sophus::SE3d::Tangent tau =
      jacobian_tangent(which, helicoid_to_sophus_tangent(take(in, "tau", 6), kSe3));
  const Sophus::SE3d::Adjoint j = is_inverse(which) ? Sophus::SE3d::leftJacobianInverse(tau)
                                                    : Sophus::SE3d::leftJacobian(tau);
  return {{"J", sophus_to_helicoid_matrix(j, kSe3)}};
}

Fields se2_out(const Sophus::SE2d& g) {
  const Eigen::Vector2d t = g.translation();
  return {{"t", {t.x(), t.y()}}, {"z", sophus_to_helicoid_so2(g.so2())}};
}

Fields se3_out(const Sophus::SE3d& g) {
  return {{"q", sophus_to_helicoid_quat(g.unit_quaternion())},
          {"x", sophus_to_helicoid_vec3(g.translation())}};
}

Sophus::SE2d se2_in(const Fields& in) { return helicoid_to_sophus_se2(take(in, "z", 2), take(in, "t", 2)); }

Sophus::SE3d se3_in(const Fields& in) { return helicoid_to_sophus_se3(take(in, "q", 4), take(in, "x", 3)); }

Sophus::Sim3d sim3_in(const Fields& in) {
  return helicoid_to_sophus_sim3(take(in, "q", 4), take(in, "x", 3), scalar(in, "sigma"));
}

// One id: its evaluation, and the fields and lengths it answers, which a refusal fills with `nan`.
struct Id {
  std::string name;
  std::function<Fields(const Fields&)> evaluate;
  std::vector<std::pair<std::string, std::size_t>> out;
};

const std::vector<Id>& table() {
  using J = Jacobian;
  static const std::vector<Id> ids = {
      {"so2_exp", [](const Fields& in) { return Fields{{"z", sophus_to_helicoid_so2(Sophus::SO2d::exp(scalar(in, "theta")))}}; },
       {{"z", 2}}},
      {"so2_log", [](const Fields& in) { return Fields{{"theta", {helicoid_to_sophus_so2(take(in, "z", 2)).log()}}}; },
       {{"theta", 1}}},
      {"se2_exp", [](const Fields& in) { return se2_out(Sophus::SE2d::exp(helicoid_to_sophus_tangent(take(in, "tau", 3), kSe2))); },
       {{"t", 2}, {"z", 2}}},
      {"se2_log", [](const Fields& in) { return Fields{{"tau", sophus_to_helicoid_tangent(se2_in(in).log(), kSe2)}}; },
       {{"tau", 3}}},
      {"se2_ad", [](const Fields& in) { return Fields{{"Ad", sophus_to_helicoid_matrix(se2_in(in).Adj(), kSe2)}}; },
       {{"Ad", 9}}},
      {"so3_exp", [](const Fields& in) { return Fields{{"q", sophus_to_helicoid_quat(Sophus::SO3d::exp(helicoid_to_sophus_vec3(take(in, "phi", 3))).unit_quaternion())}}; },
       {{"q", 4}}},
      {"so3_log", [](const Fields& in) { return Fields{{"phi", sophus_to_helicoid_vec3(helicoid_to_sophus_so3(take(in, "q", 4)).log())}}; },
       {{"phi", 3}}},
      {"so3_act", [](const Fields& in) { return Fields{{"Rp", sophus_to_helicoid_vec3(helicoid_to_sophus_so3(take(in, "q", 4)) * helicoid_to_sophus_vec3(take(in, "p", 3)))}}; },
       {{"Rp", 3}}},
      // `SO3(R)`: Eigen's quaternion of the matrix, after Sophus's orthogonality check.
      {"so3_from_matrix", [](const Fields& in) {
         const Eigen::Matrix3d r = Eigen::Map<const Eigen::Matrix3d>(take(in, "R", 9).data());
         return Fields{{"q", sophus_to_helicoid_quat(Sophus::SO3d(r).unit_quaternion())}};
       },
       {{"q", 4}}},
      {"so3_jl", [](const Fields& in) { return so3_jacobian(in, J::Left); }, {{"J", 9}}},
      {"so3_jr", [](const Fields& in) { return so3_jacobian(in, J::Right); }, {{"J", 9}}},
      {"so3_jl_inv", [](const Fields& in) { return so3_jacobian(in, J::LeftInverse); }, {{"J", 9}}},
      {"so3_jr_inv", [](const Fields& in) { return so3_jacobian(in, J::RightInverse); }, {{"J", 9}}},
      // `γ(X₀, X₁, t) = X₀ Exp(t Log(X₀⁻¹ X₁))` (docs/NUMERICS.md §10) is `Sophus::interpolate`.
      {"so3_geodesic", [](const Fields& in) {
         const auto q = Sophus::interpolate(helicoid_to_sophus_so3(take(in, "q0", 4)),
                                            helicoid_to_sophus_so3(take(in, "q1", 4)), scalar(in, "t"));
         return Fields{{"q", sophus_to_helicoid_quat(q.unit_quaternion())}};
       },
       {{"q", 4}}},
      {"se3_geodesic", [](const Fields& in) {
         return se3_out(Sophus::interpolate(helicoid_to_sophus_se3(take(in, "q0", 4), take(in, "x0", 3)),
                                            helicoid_to_sophus_se3(take(in, "q1", 4), take(in, "x1", 3)),
                                            scalar(in, "t")));
       },
       {{"q", 4}, {"x", 3}}},
      {"sen3_exp_n1", [](const Fields& in) { return se3_out(Sophus::SE3d::exp(helicoid_to_sophus_tangent(take(in, "tau", 6), kSe3))); },
       {{"q", 4}, {"x", 3}}},
      {"sen3_log_n1", [](const Fields& in) { return Fields{{"tau", sophus_to_helicoid_tangent(se3_in(in).log(), kSe3)}}; },
       {{"tau", 6}}},
      {"sen3_ad_n1", [](const Fields& in) { return Fields{{"Ad", sophus_to_helicoid_matrix(se3_in(in).Adj(), kSe3)}}; },
       {{"Ad", 36}}},
      {"sen3_jl_n1", [](const Fields& in) { return se3_jacobian(in, J::Left); }, {{"J", 36}}},
      {"sen3_jr_n1", [](const Fields& in) { return se3_jacobian(in, J::Right); }, {{"J", 36}}},
      {"sen3_jl_inv_n1", [](const Fields& in) { return se3_jacobian(in, J::LeftInverse); }, {{"J", 36}}},
      {"sen3_jr_inv_n1", [](const Fields& in) { return se3_jacobian(in, J::RightInverse); }, {{"J", 36}}},
      {"sim3_exp", [](const Fields& in) {
         const auto x = sophus_to_helicoid_sim3(Sophus::Sim3d::exp(helicoid_to_sophus_tangent(take(in, "tau", 7), kSim3)));
         return Fields{{"q", x.q}, {"x", x.x}, {"sigma", {x.sigma}}};
       },
       {{"q", 4}, {"x", 3}, {"sigma", 1}}},
      {"sim3_log", [](const Fields& in) { return Fields{{"tau", sophus_to_helicoid_tangent(sim3_in(in).log(), kSim3)}}; },
       {{"tau", 7}}},
      {"sim3_ad", [](const Fields& in) { return Fields{{"Ad", sophus_to_helicoid_matrix(sim3_in(in).Adj(), kSim3)}}; },
       {{"Ad", 49}}},
      {"sim3_act", [](const Fields& in) { return Fields{{"Xp", sophus_to_helicoid_vec3(sim3_in(in) * helicoid_to_sophus_vec3(take(in, "p", 3)))}}; },
       {{"Xp", 3}}},
  };
  return ids;
}

const Id* find(const std::string& name) {
  for (const Id& id : table()) {
    if (id.name == name) return &id;
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

// The answer to one record: Sophus's, its shape checked, or `nan` throughout if Sophus refuses.
Fields answer(const Id& id, const Fields& in, std::size_t& refused) {
  Fields out;
  try {
    out = id.evaluate(in);
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

// Answers every line of `path` into `target`; returns the number of refused records.
std::size_t answer_file(const Id& id, const std::string& path, const std::string& target) {
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
      const std::uint64_t record_id = record.at("id").get<std::uint64_t>();
      out << line_of(record_id, answer(id, inputs(record.at("in")), refused));
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

int run(const std::vector<std::string>& args) {
  if (args.size() == 1 && args[0] == "--version") {
    std::cout << pin() << '\n';
    return std::cout.flush() ? 0 : 1;
  }
  if (args.size() < 3 || args[0] != "--out") {
    std::cerr << "usage: sophus_runner --version | --out DIR FILE.jsonl…\n";
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
    const Id* id = find(stem);
    if (id == nullptr) {
      ++skipped;
      continue;
    }
    refused += answer_file(*id, path, dir + "/" + stem + ".jsonl");
    answered.push_back(stem);
  }
  std::cerr << "answered " << answered.size() << " ids; " << refused << " records refused (nan); " << skipped
            << " files of other ids skipped\n";
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
