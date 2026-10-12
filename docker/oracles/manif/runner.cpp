// The manif oracle (docs/PHASE1.md §7, docs/decisions/0075): reads corpus JSONL and writes one answer
// file per function id it answers. It never reads a reference and never scores;
// `cargo xtask conformance --oracle manif` does, with the harness's exact metric.
//
//   manif_runner --version               the pin, one line: the result rows' subject_version
//   manif_runner --out DIR FILE.jsonl…   DIR/<id>.jsonl for each FILE whose id is answered
//
// An answer line is `{"id":N,"out":{"<field>":["<hex float>",…]}}`, fields flat in the reference's
// order (matrices column-major). A record manif refuses (it throws `manif::invalid_argument` or
// `manif::runtime_error`) is answered with `nan` in every field: answered, so no stratum is the
// runner's choice, and non-finite, so it never sets a bar.
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

#include "convert.hpp"
#include "hexfloat.hpp"

#ifndef ORACLE_MANIF
#error "ORACLE_MANIF (the pinned manif version and commit) is set by the Dockerfile"
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

// Which Jacobian an id asks for: manif has all four, each its own formula (`ljac`, `rjac`,
// `ljacinv`, `rjacinv`), so none is answered through another.
enum class Jacobian { Left, Right, LeftInverse, RightInverse };

template <class T>
auto jacobian(const T& tau, Jacobian which) -> decltype(tau.ljac()) {
  switch (which) {
    case Jacobian::Left:
      return tau.ljac();
    case Jacobian::Right:
      return tau.rjac();
    case Jacobian::LeftInverse:
      return tau.ljacinv();
    case Jacobian::RightInverse:
      return tau.rjacinv();
  }
  throw std::logic_error("no such Jacobian");
}

Fields so3_jacobian(const Fields& in, Jacobian which) {
  const manif::SO3Tangentd phi(helicoid_to_manif_vec3(take(in, "phi", 3)));
  return {{"J", manif_to_helicoid_matrix3(jacobian(phi, which))}};
}

Fields se2_jacobian(const Fields& in, Jacobian which) {
  const manif::SE2Tangentd tau(helicoid_to_manif_tangent(take(in, "tau", 3), kSe2));
  return {{"J", manif_to_helicoid_matrix(jacobian(tau, which), kSe2)}};
}

Fields se3_jacobian(const Fields& in, Jacobian which) {
  const manif::SE3Tangentd tau(helicoid_to_manif_tangent(take(in, "tau", 6), kSe3));
  return {{"J", manif_to_helicoid_matrix(jacobian(tau, which), kSe3)}};
}

Fields se23_jacobian(const Fields& in, Jacobian which) {
  const manif::SE_2_3Tangentd tau(helicoid_to_manif_tangent(take(in, "tau", 9), kSe23));
  return {{"J", manif_to_helicoid_matrix(jacobian(tau, which), kSe23)}};
}

Fields se2_out(const manif::SE2d& g) {
  const HelicoidSe2 x = manif_to_helicoid_se2(g);
  return {{"t", x.t}, {"z", x.z}};
}

Fields sen3_out(const HelicoidSen3& x) { return {{"q", x.q}, {"x", x.x}}; }

// A chart's `local` `[φ; ρ]`, held as two fields.
Fields local_out(const Vec& tau) { return {{"phi", {tau[0], tau[1], tau[2]}}, {"rho", {tau[3], tau[4], tau[5]}}}; }

manif::SE2d se2_in(const Fields& in) { return helicoid_to_manif_se2(take(in, "z", 2), take(in, "t", 2)); }

manif::SE3d se3_in(const Fields& in, const std::string& q, const std::string& x) {
  return helicoid_to_manif_se3(take(in, q, 4), take(in, x, 3));
}

manif::SE_2_3d se23_in(const Fields& in) { return helicoid_to_manif_se23(take(in, "q", 4), take(in, "x", 6)); }

SO3xR3 so3xr3_in(const Fields& in, const std::string& q, const std::string& x) {
  return helicoid_to_manif_so3xr3(take(in, q, 4), take(in, x, 3));
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
      {"so2_exp", [](const Fields& in) { return Fields{{"z", manif_to_helicoid_so2(manif::SO2Tangentd(scalar(in, "theta")).exp())}}; },
       {{"z", 2}}},
      {"so2_log", [](const Fields& in) { return Fields{{"theta", {helicoid_to_manif_so2(take(in, "z", 2)).log().angle()}}}; },
       {{"theta", 1}}},
      {"se2_exp", [](const Fields& in) { return se2_out(manif::SE2Tangentd(helicoid_to_manif_tangent(take(in, "tau", 3), kSe2)).exp()); },
       {{"t", 2}, {"z", 2}}},
      {"se2_log", [](const Fields& in) { return Fields{{"tau", manif_to_helicoid_tangent(se2_in(in).log().coeffs(), kSe2)}}; },
       {{"tau", 3}}},
      {"se2_ad", [](const Fields& in) { return Fields{{"Ad", manif_to_helicoid_matrix(se2_in(in).adj(), kSe2)}}; },
       {{"Ad", 9}}},
      {"se2_jl", [](const Fields& in) { return se2_jacobian(in, J::Left); }, {{"J", 9}}},
      {"se2_jr", [](const Fields& in) { return se2_jacobian(in, J::Right); }, {{"J", 9}}},
      {"se2_jl_inv", [](const Fields& in) { return se2_jacobian(in, J::LeftInverse); }, {{"J", 9}}},
      {"se2_jr_inv", [](const Fields& in) { return se2_jacobian(in, J::RightInverse); }, {{"J", 9}}},
      {"so3_exp", [](const Fields& in) { return Fields{{"q", manif_to_helicoid_quat(manif::SO3Tangentd(helicoid_to_manif_vec3(take(in, "phi", 3))).exp().quat())}}; },
       {{"q", 4}}},
      {"so3_log", [](const Fields& in) { return Fields{{"phi", manif_to_helicoid_vec3(helicoid_to_manif_so3(take(in, "q", 4)).log().coeffs())}}; },
       {{"phi", 3}}},
      {"so3_act", [](const Fields& in) { return Fields{{"Rp", manif_to_helicoid_vec3(helicoid_to_manif_so3(take(in, "q", 4)).act(helicoid_to_manif_vec3(take(in, "p", 3))))}}; },
       {{"Rp", 3}}},
      {"so3_jl", [](const Fields& in) { return so3_jacobian(in, J::Left); }, {{"J", 9}}},
      {"so3_jr", [](const Fields& in) { return so3_jacobian(in, J::Right); }, {{"J", 9}}},
      {"so3_jl_inv", [](const Fields& in) { return so3_jacobian(in, J::LeftInverse); }, {{"J", 9}}},
      {"so3_jr_inv", [](const Fields& in) { return so3_jacobian(in, J::RightInverse); }, {{"J", 9}}},
      // `γ(X₀, X₁, t) = X₀ Exp(t Log(X₀⁻¹ X₁))` (docs/NUMERICS.md §10) is `interpolate_slerp`,
      // `X₀.rplus(X₁.rminus(X₀) · t)`.
      {"so3_geodesic", [](const Fields& in) {
         const manif::SO3d q = manif::interpolate_slerp(helicoid_to_manif_so3(take(in, "q0", 4)),
                                                        helicoid_to_manif_so3(take(in, "q1", 4)), scalar(in, "t"));
         return Fields{{"q", manif_to_helicoid_quat(q.quat())}};
       },
       {{"q", 4}}},
      {"se3_geodesic", [](const Fields& in) {
         return sen3_out(manif_to_helicoid_se3(
             manif::interpolate_slerp(se3_in(in, "q0", "x0"), se3_in(in, "q1", "x1"), scalar(in, "t"))));
       },
       {{"q", 4}, {"x", 3}}},
      // The screw chart is SE(3)'s `rplus` and `rminus`: `X Exp(τ)`, `Log(X⁻¹ Y)`.
      {"se3_screw_retract", [](const Fields& in) {
         const manif::SE3Tangentd tau(helicoid_to_manif_tangent(take(in, "tau", 6), kSe3));
         return sen3_out(manif_to_helicoid_se3(se3_in(in, "q0", "x0").rplus(tau)));
       },
       {{"q", 4}, {"x", 3}}},
      {"se3_screw_local", [](const Fields& in) {
         return local_out(manif_to_helicoid_tangent(se3_in(in, "q1", "x1").rminus(se3_in(in, "q0", "x0")).coeffs(), kSe3));
       },
       {{"phi", 3}, {"rho", 3}}},
      // The world chart is `Bundle<SO3, R3>`'s `rplus` and `rminus`: `(q Exp φ, x + ρ)`,
      // `(Log(q₀* q₁), x₁ − x₀)`, its tangent `[φ; ρ]` already helicoid's order.
      {"se3_world_retract", [](const Fields& in) {
         const SO3xR3::Tangent tau(Eigen::Map<const SO3xR3::Tangent::DataType>(take(in, "tau", 6).data()));
         return sen3_out(manif_to_helicoid_so3xr3(so3xr3_in(in, "q0", "x0").rplus(tau)));
       },
       {{"q", 4}, {"x", 3}}},
      {"se3_world_local", [](const Fields& in) {
         const SO3xR3::Tangent::DataType tau = so3xr3_in(in, "q1", "x1").rminus(so3xr3_in(in, "q0", "x0")).coeffs();
         return local_out(Vec(tau.data(), tau.data() + 6));
       },
       {{"phi", 3}, {"rho", 3}}},
      {"sen3_exp_n1", [](const Fields& in) { return sen3_out(manif_to_helicoid_se3(manif::SE3Tangentd(helicoid_to_manif_tangent(take(in, "tau", 6), kSe3)).exp())); },
       {{"q", 4}, {"x", 3}}},
      {"sen3_log_n1", [](const Fields& in) { return Fields{{"tau", manif_to_helicoid_tangent(se3_in(in, "q", "x").log().coeffs(), kSe3)}}; },
       {{"tau", 6}}},
      {"sen3_ad_n1", [](const Fields& in) { return Fields{{"Ad", manif_to_helicoid_matrix(se3_in(in, "q", "x").adj(), kSe3)}}; },
       {{"Ad", 36}}},
      {"sen3_jl_n1", [](const Fields& in) { return se3_jacobian(in, J::Left); }, {{"J", 36}}},
      {"sen3_jr_n1", [](const Fields& in) { return se3_jacobian(in, J::Right); }, {{"J", 36}}},
      {"sen3_jl_inv_n1", [](const Fields& in) { return se3_jacobian(in, J::LeftInverse); }, {{"J", 36}}},
      {"sen3_jr_inv_n1", [](const Fields& in) { return se3_jacobian(in, J::RightInverse); }, {{"J", 36}}},
      {"sen3_exp_n2", [](const Fields& in) { return sen3_out(manif_to_helicoid_se23(manif::SE_2_3Tangentd(helicoid_to_manif_tangent(take(in, "tau", 9), kSe23)).exp())); },
       {{"q", 4}, {"x", 6}}},
      {"sen3_log_n2", [](const Fields& in) { return Fields{{"tau", manif_to_helicoid_tangent(se23_in(in).log().coeffs(), kSe23)}}; },
       {{"tau", 9}}},
      {"sen3_ad_n2", [](const Fields& in) { return Fields{{"Ad", manif_to_helicoid_matrix(se23_in(in).adj(), kSe23)}}; },
       {{"Ad", 81}}},
      {"sen3_jl_n2", [](const Fields& in) { return se23_jacobian(in, J::Left); }, {{"J", 81}}},
      {"sen3_jr_n2", [](const Fields& in) { return se23_jacobian(in, J::Right); }, {{"J", 81}}},
      {"sen3_jl_inv_n2", [](const Fields& in) { return se23_jacobian(in, J::LeftInverse); }, {{"J", 81}}},
      {"sen3_jr_inv_n2", [](const Fields& in) { return se23_jacobian(in, J::RightInverse); }, {{"J", 81}}},
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

// The answer to one record: manif's, its shape checked, or `nan` throughout if manif refuses.
Fields answer(const Id& id, const Fields& in, std::size_t& refused) {
  Fields out;
  try {
    out = id.evaluate(in);
  } catch (const manif::invalid_argument&) {
    ++refused;
    for (const auto& [field, len] : id.out) out[field] = Vec(len, std::nan(""));
  } catch (const manif::runtime_error&) {
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
  text << "manif@" << ORACLE_MANIF << " eigen@" << EIGEN_WORLD_VERSION << '.' << EIGEN_MAJOR_VERSION << '.'
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
    std::cerr << "usage: manif_runner --version | --out DIR FILE.jsonl…\n";
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
    std::cerr << "manif_runner: " << e.what() << '\n';
    return 1;
  }
}
