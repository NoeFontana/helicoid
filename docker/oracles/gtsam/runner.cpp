// The GTSAM C++ oracle (docs/PHASE1.md §7, docs/decisions/0074): reads corpus JSONL and writes one
// answer file per function id it answers. It never reads a reference and never scores;
// `cargo xtask conformance --oracle gtsam` does, with the harness's exact metric.
//
//   gtsam_runner --version               the pin, one line: the result rows' subject_version
//   gtsam_runner --out DIR FILE.jsonl…   DIR/<id>.jsonl for each FILE whose id is answered
//
// An answer line is `{"id":N,"out":{"<field>":["<hex float>",…]}}`, fields flat in the reference's
// order (matrices column-major). GTSAM refuses no record of the corpus, so an exception is a failed
// run, not an answer.
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

#include <gtsam/config.h>
#include <gtsam/geometry/Kernel.h>
#include <gtsam/geometry/SO3.h>
#include <nlohmann/json.hpp>

#include "convert.hpp"
#include "hexfloat.hpp"

#ifndef ORACLE_GTSAM
#error "ORACLE_GTSAM (the pinned GTSAM version and commit) is set by the Dockerfile"
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

// Which Jacobian an id asks for. GTSAM's `ExpmapDerivative` is the right Jacobian and its
// `LogmapDerivative` of a tangent the right one's inverse; the left ones are answered at the negated
// tangent, `J_l(τ) = J_r(−τ)` (docs/NUMERICS.md §1).
enum class Jacobian { Left, Right, LeftInverse, RightInverse };

template <class T>
T jacobian_tangent(Jacobian which, const T& tau) {
  return which == Jacobian::Left || which == Jacobian::LeftInverse ? T(-tau) : tau;
}

bool is_inverse(Jacobian which) {
  return which == Jacobian::LeftInverse || which == Jacobian::RightInverse;
}

// `G` is a GTSAM class with static `ExpmapDerivative` and `LogmapDerivative` of a tangent.
template <class G, class T>
auto right_jacobian(Jacobian which, const T& tau) {
  const T at = jacobian_tangent(which, tau);
  return is_inverse(which) ? G::LogmapDerivative(at) : G::ExpmapDerivative(at);
}

Fields so3_jacobian(const Fields& in, Jacobian which) {
  const Eigen::Vector3d phi = helicoid_to_gtsam_vec3(take(in, "phi", 3));
  return {{"J", gtsam_to_helicoid_matrix(right_jacobian<gtsam::Rot3>(which, phi))}};
}

Fields se2_jacobian(const Fields& in, Jacobian which) {
  const Eigen::Vector3d tau = helicoid_to_gtsam_se2_tangent(take(in, "tau", 3));
  return {{"J", gtsam_to_helicoid_se2_matrix(right_jacobian<gtsam::Pose2>(which, tau))}};
}

// SE_N(3): `Pose3` at N = 1, `ExtendedPose3<N>` above it.
template <int K>
struct SeN3;

template <>
struct SeN3<1> {
  using G = gtsam::Pose3;
  static G in(const Fields& in) { return helicoid_to_gtsam_pose3(take(in, "q", 4), take(in, "x", 3)); }
  static Fields out(const G& g) {
    return {{"q", gtsam_to_helicoid_quat(g.rotation())}, {"x", gtsam_to_helicoid_vec3(g.translation())}};
  }
};

template <int K>
struct SeN3 {
  using G = gtsam::ExtendedPose3<K>;
  static G in(const Fields& in) {
    return helicoid_to_gtsam_extended_pose3<K>(take(in, "q", 4), take(in, "x", 3 * K));
  }
  static Fields out(const G& g) {
    return {{"q", gtsam_to_helicoid_quat(g.rotation())}, {"x", gtsam_to_helicoid_matrix(g.xMatrix())}};
  }
};

template <int K>
Fields sen3_exp(const Fields& in) {
  using G = typename SeN3<K>::G;
  return SeN3<K>::out(G::Expmap(helicoid_to_gtsam_tangent<3 + 3 * K>(take(in, "tau", 3 + 3 * K))));
}

template <int K>
Fields sen3_log(const Fields& in) {
  using G = typename SeN3<K>::G;
  return {{"tau", gtsam_to_helicoid_tangent<3 + 3 * K>(G::Logmap(SeN3<K>::in(in)))}};
}

template <int K>
Fields sen3_ad(const Fields& in) {
  return {{"Ad", gtsam_to_helicoid_matrix(SeN3<K>::in(in).AdjointMap())}};
}

template <int K>
Fields sen3_jacobian(const Fields& in, Jacobian which) {
  using G = typename SeN3<K>::G;
  const Eigen::Matrix<double, 3 + 3 * K, 1> tau = helicoid_to_gtsam_tangent<3 + 3 * K>(take(in, "tau", 3 + 3 * K));
  const Eigen::Matrix<double, 3 + 3 * K, 3 + 3 * K> j = right_jacobian<G>(which, tau);
  return {{"J", gtsam_to_helicoid_matrix(j)}};
}

Fields se2_out(const gtsam::Pose2& g) {
  return {{"t", {g.t().x(), g.t().y()}}, {"z", gtsam_to_helicoid_rot2(g.r())}};
}

gtsam::Pose2 se2_in(const Fields& in) { return helicoid_to_gtsam_pose2(take(in, "z", 2), take(in, "t", 2)); }

gtsam::Similarity3 sim3_in(const Fields& in) {
  return helicoid_to_gtsam_similarity3(take(in, "q", 4), take(in, "x", 3), scalar(in, "sigma"));
}

// One id: its evaluation, and the fields and lengths it answers.
struct Id {
  std::string name;
  std::function<Fields(const Fields&)> evaluate;
  std::vector<std::pair<std::string, std::size_t>> out;
};

template <int K>
void add_sen3(std::vector<Id>& ids) {
  using J = Jacobian;
  const std::string n = "_n" + std::to_string(K);
  const std::size_t d = 3 + 3 * K;
  ids.push_back({"sen3_exp" + n, sen3_exp<K>, {{"q", 4}, {"x", 3 * K}}});
  ids.push_back({"sen3_log" + n, sen3_log<K>, {{"tau", d}}});
  ids.push_back({"sen3_ad" + n, sen3_ad<K>, {{"Ad", d * d}}});
  ids.push_back({"sen3_jl" + n, [](const Fields& in) { return sen3_jacobian<K>(in, J::Left); }, {{"J", d * d}}});
  ids.push_back({"sen3_jr" + n, [](const Fields& in) { return sen3_jacobian<K>(in, J::Right); }, {{"J", d * d}}});
  ids.push_back({"sen3_jl_inv" + n, [](const Fields& in) { return sen3_jacobian<K>(in, J::LeftInverse); }, {{"J", d * d}}});
  ids.push_back({"sen3_jr_inv" + n, [](const Fields& in) { return sen3_jacobian<K>(in, J::RightInverse); }, {{"J", d * d}}});
}

std::vector<Id> make_table() {
  using J = Jacobian;
  std::vector<Id> ids = {
      {"so2_exp", [](const Fields& in) { return Fields{{"z", gtsam_to_helicoid_rot2(gtsam::Rot2::Expmap(gtsam::Vector1(scalar(in, "theta"))))}}; },
       {{"z", 2}}},
      {"so2_log", [](const Fields& in) { return Fields{{"theta", {gtsam::Rot2::Logmap(helicoid_to_gtsam_rot2(take(in, "z", 2)))(0)}}}; },
       {{"theta", 1}}},
      {"se2_exp", [](const Fields& in) { return se2_out(gtsam::Pose2::Expmap(helicoid_to_gtsam_se2_tangent(take(in, "tau", 3)))); },
       {{"t", 2}, {"z", 2}}},
      {"se2_log", [](const Fields& in) { return Fields{{"tau", gtsam_to_helicoid_se2_tangent(gtsam::Pose2::Logmap(se2_in(in)))}}; },
       {{"tau", 3}}},
      {"se2_ad", [](const Fields& in) { return Fields{{"Ad", gtsam_to_helicoid_se2_matrix(se2_in(in).AdjointMap())}}; },
       {{"Ad", 9}}},
      {"se2_jl", [](const Fields& in) { return se2_jacobian(in, J::Left); }, {{"J", 9}}},
      {"se2_jr", [](const Fields& in) { return se2_jacobian(in, J::Right); }, {{"J", 9}}},
      {"se2_jl_inv", [](const Fields& in) { return se2_jacobian(in, J::LeftInverse); }, {{"J", 9}}},
      {"se2_jr_inv", [](const Fields& in) { return se2_jacobian(in, J::RightInverse); }, {{"J", 9}}},
      {"so3_exp", [](const Fields& in) { return Fields{{"q", gtsam_to_helicoid_quat(gtsam::Rot3::Expmap(helicoid_to_gtsam_vec3(take(in, "phi", 3))))}}; },
       {{"q", 4}}},
      {"so3_log", [](const Fields& in) { return Fields{{"phi", gtsam_to_helicoid_vec3(gtsam::Rot3::Logmap(helicoid_to_gtsam_rot3(take(in, "q", 4))))}}; },
       {{"phi", 3}}},
      {"so3_act", [](const Fields& in) { return Fields{{"Rp", gtsam_to_helicoid_vec3(helicoid_to_gtsam_rot3(take(in, "q", 4)).rotate(helicoid_to_gtsam_vec3(take(in, "p", 3))))}}; },
       {{"Rp", 3}}},
      // `Rot3(R).toQuaternion()`: `Rot3` keeps the matrix as given, so this is Eigen's quaternion
      // of it.
      {"so3_from_matrix", [](const Fields& in) { return Fields{{"q", gtsam_to_helicoid_quat(gtsam::Rot3(helicoid_to_gtsam_matrix3(take(in, "R", 9))))}}; },
       {{"q", 4}}},
      {"so3_jl", [](const Fields& in) { return so3_jacobian(in, J::Left); }, {{"J", 9}}},
      {"so3_jr", [](const Fields& in) { return so3_jacobian(in, J::Right); }, {{"J", 9}}},
      {"so3_jl_inv", [](const Fields& in) { return so3_jacobian(in, J::LeftInverse); }, {{"J", 9}}},
      {"so3_jr_inv", [](const Fields& in) { return so3_jacobian(in, J::RightInverse); }, {{"J", 9}}},
      // `Γ₂(φ) = Σ Wⁿ/(n+2)!` (docs/NUMERICS.md §7) is the left `Gamma` kernel of GTSAM's
      // `DexpFunctor`, `½I + C Ω + E Ω²`.
      {"so3_gamma2", [](const Fields& in) {
         const gtsam::so3::DexpFunctor local(helicoid_to_gtsam_vec3(take(in, "phi", 3)));
         return Fields{{"G", gtsam_to_helicoid_matrix(local.Gamma().left())}};
       },
       {{"G", 9}}},
      // `γ(X₀, X₁, t) = X₀ Exp(t Log(X₀⁻¹ X₁))` (docs/NUMERICS.md §10) is GTSAM's `slerp`, its
      // `interpolate` through `Expmap` and `Logmap`.
      {"so3_geodesic", [](const Fields& in) {
         const gtsam::Rot3 r = helicoid_to_gtsam_rot3(take(in, "q0", 4)).slerp(scalar(in, "t"), helicoid_to_gtsam_rot3(take(in, "q1", 4)));
         return Fields{{"q", gtsam_to_helicoid_quat(r)}};
       },
       {{"q", 4}}},
      {"se3_geodesic", [](const Fields& in) {
         const gtsam::Pose3 x0 = helicoid_to_gtsam_pose3(take(in, "q0", 4), take(in, "x0", 3));
         const gtsam::Pose3 x1 = helicoid_to_gtsam_pose3(take(in, "q1", 4), take(in, "x1", 3));
         return SeN3<1>::out(x0.slerp(scalar(in, "t"), x1));
       },
       {{"q", 4}, {"x", 3}}},
      // `Pose3`'s chart is `Expmap` (`GTSAM_POSE3_EXPMAP`, the default): `retract` is `X Exp(τ)` and
      // `localCoordinates` is `Log(X⁻¹ Y)`, the screw chart.
      {"se3_screw_retract", [](const Fields& in) {
         const gtsam::Pose3 x0 = helicoid_to_gtsam_pose3(take(in, "q0", 4), take(in, "x0", 3));
         return SeN3<1>::out(x0.retract(helicoid_to_gtsam_tangent<6>(take(in, "tau", 6))));
       },
       {{"q", 4}, {"x", 3}}},
      {"se3_screw_local", [](const Fields& in) {
         const gtsam::Pose3 x0 = helicoid_to_gtsam_pose3(take(in, "q0", 4), take(in, "x0", 3));
         const gtsam::Pose3 x1 = helicoid_to_gtsam_pose3(take(in, "q1", 4), take(in, "x1", 3));
         const gtsam::Vector6 tau = x0.localCoordinates(x1);
         return Fields{{"phi", {tau(0), tau(1), tau(2)}}, {"rho", {tau(3), tau(4), tau(5)}}};
       },
       {{"phi", 3}, {"rho", 3}}},
      {"sim3_exp", [](const Fields& in) {
         const auto x = gtsam_to_helicoid_sim3(gtsam::Similarity3::Expmap(helicoid_to_gtsam_tangent<7>(take(in, "tau", 7))));
         return Fields{{"q", x.q}, {"x", x.x}, {"sigma", {x.sigma}}};
       },
       {{"q", 4}, {"x", 3}, {"sigma", 1}}},
      {"sim3_log", [](const Fields& in) { return Fields{{"tau", gtsam_to_helicoid_tangent<7>(gtsam::Similarity3::Logmap(sim3_in(in)))}}; },
       {{"tau", 7}}},
      {"sim3_ad", [](const Fields& in) { return Fields{{"Ad", gtsam_to_helicoid_matrix(sim3_in(in).AdjointMap())}}; },
       {{"Ad", 49}}},
      {"sim3_act", [](const Fields& in) { return Fields{{"Xp", gtsam_to_helicoid_vec3(sim3_in(in).transformFrom(helicoid_to_gtsam_vec3(take(in, "p", 3))))}}; },
       {{"Xp", 3}}},
  };
  add_sen3<1>(ids);
  add_sen3<2>(ids);
  add_sen3<3>(ids);
  return ids;
}

const std::vector<Id>& table() {
  static const std::vector<Id> ids = make_table();
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

// The answer to one record, its shape checked.
Fields answer(const Id& id, const Fields& in) {
  const Fields out = id.evaluate(in);
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

// Answers every line of `path` into `target`.
void answer_file(const Id& id, const std::string& path, const std::string& target) {
  std::ifstream file(path);
  if (!file) throw std::runtime_error(path + ": cannot read");
  std::ofstream out(target);
  if (!out) throw std::runtime_error(target + ": cannot write");
  std::string text;
  std::size_t number = 0;
  while (std::getline(file, text)) {
    ++number;
    try {
      const auto record = nlohmann::json::parse(text);
      const std::uint64_t record_id = record.at("id").get<std::uint64_t>();
      out << line_of(record_id, answer(id, inputs(record.at("in"))));
    } catch (const std::exception& e) {
      throw std::runtime_error(path + ": line " + std::to_string(number) + ": " + e.what());
    }
  }
  if (!out.flush()) throw std::runtime_error(target + ": write failed");
}

std::string pin() {
  std::ostringstream text;
  text << "gtsam@" << ORACLE_GTSAM << " eigen@" << EIGEN_WORLD_VERSION << '.' << EIGEN_MAJOR_VERSION << '.'
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
    std::cerr << "usage: gtsam_runner --version | --out DIR FILE.jsonl…\n";
    return 2;
  }
  const std::string dir = args[1];
  std::size_t answered = 0, skipped = 0;
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
    answer_file(*id, path, dir + "/" + stem + ".jsonl");
    ++answered;
  }
  std::cerr << "answered " << answered << " ids; " << skipped << " files of other ids skipped\n";
  return 0;
}

}  // namespace
}  // namespace oracle

int main(int argc, char** argv) {
  try {
    return oracle::run(std::vector<std::string>(argv + 1, argv + argc));
  } catch (const std::exception& e) {
    std::cerr << "gtsam_runner: " << e.what() << '\n';
    return 1;
  }
}
