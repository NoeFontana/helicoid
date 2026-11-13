// The basalt-headers oracle (docs/PHASE1.md §7, docs/decisions/0116 decision 4): reads corpus JSONL
// and writes one answer file per spline id it answers. It never reads a reference and never scores;
// `cargo xtask conformance --oracle basalt` does, with the harness's exact metric.
//
//   basalt_runner --version                 the pin, one line: the result rows' subject_version
//   basalt_runner --out DIR FILE.jsonl…     DIR/<id>.jsonl for each FILE whose id is answered
//
// Binary64 only: `CeresSplineHelper`'s weights are `double` whatever its knots' scalar, and
// `So3Spline`'s Jacobians call `Sophus::SO3d`.
//
// Each id is one function of `spline.hpp`. A record off the integer-nanosecond grid, or one basalt
// refuses, is answered `nan`: answered, so no stratum is the runner's choice, and non-finite, so it
// never sets a bar (docs/decisions/0072 decision 5).
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
#include "spline.hpp"

#ifndef ORACLE_BASALT
#error "ORACLE_BASALT (the pinned basalt-headers commit) is set by the Dockerfile"
#endif
#ifndef ORACLE_SOPHUS
#error "ORACLE_SOPHUS (the pinned Sophus version and commit) is set by the Dockerfile"
#endif

namespace oracle {
namespace {

// An answer's fields and their lengths, which a refusal fills with `nan`.
using Shape = std::vector<std::pair<std::string, std::size_t>>;

// --- the ids -------------------------------------------------------------------------------------

struct Id {
  std::string name;
  std::function<Fields(const Fields&)> evaluate;
  // The answer's shape for this record, from its order.
  std::function<Shape(const Fields&)> shape;
};

const std::vector<Id>& table() {
  const auto fixed = [](Shape s) { return [s](const Fields&) { return s; }; };
  const auto per_control = [](const char* key, std::size_t per, const char* field, std::size_t each) {
    return [=](const Fields& in) { return Shape{{field, order_of(in, key, per) * each}}; };
  };
  static const std::vector<Id> ids = {
      {"so3_spline_value", so3_rate<0>, fixed({{"q", 4}})},
      {"so3_spline_vel", so3_rate<1>, fixed({{"omega", 3}})},
      {"so3_spline_acc", so3_rate<2>, fixed({{"alpha", 3}})},
      {"so3_spline_value_jac", so3_jacobian<0>, per_control("q", 4, "jac", 9)},
      {"so3_spline_vel_jac", so3_jacobian<1>, per_control("q", 4, "jac", 9)},
      {"so3_spline_acc_jac", so3_jacobian<2>, per_control("q", 4, "jac", 9)},
      {"rn_spline_value", rn_rate<0>, fixed({{"p", 3}})},
      {"rn_spline_vel", rn_rate<1>, fixed({{"v", 3}})},
      {"rn_spline_acc", rn_rate<2>, fixed({{"a", 3}})},
      {"rn_spline_value_jac", rn_jacobian<0>, per_control("p", 3, "k", 1)},
      {"rn_spline_vel_jac", rn_jacobian<1>, per_control("p", 3, "k", 1)},
      {"rn_spline_acc_jac", rn_jacobian<2>, per_control("p", 3, "k", 1)},
      {"se3_spline_value", se3_rate<0>, fixed({{"q", 4}, {"x", 3}})},
      {"se3_spline_vel", se3_rate<1>, fixed({{"twist", 6}})},
      {"se3_spline_acc", se3_rate<2>, fixed({{"twist", 6}})},
      {"rn_spline_gram", rn_gram,
       [](const Fields& in) {
         const auto k = static_cast<std::size_t>(scalar(in, "k"));
         return Shape{{"g", k * k}};
       }},
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

// A record's `in` object: each entry a hex-float string or an array of them.
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

struct Counts {
  std::size_t off_grid = 0, refused = 0;
};

// The answer to one record: basalt's, its shape checked, or `nan` throughout where it is refused.
Fields answer(const Id& id, const Fields& in, Counts& counts) {
  const Shape shape = id.shape(in);
  Fields out;
  try {
    out = id.evaluate(in);
  } catch (const Refused& e) {
    ++(dynamic_cast<const OffGrid*>(&e) ? counts.off_grid : counts.refused);
    for (const auto& [field, len] : shape) out[field] = Vec(len, std::nan(""));
  }
  if (out.size() != shape.size()) throw std::logic_error(id.name + ": wrong output fields");
  for (const auto& [field, len] : shape) {
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
void answer_file(const Id& id, const std::string& path, const std::string& target, Counts& counts) {
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
      out << line_of(record_id, answer(id, inputs(record.at("in")), counts));
    } catch (const std::exception& e) {
      throw std::runtime_error(path + ": line " + std::to_string(number) + ": " + e.what());
    }
  }
  if (!out.flush()) throw std::runtime_error(target + ": write failed");
}

std::string pin() {
  std::ostringstream text;
  text << "basalt-headers@" << ORACLE_BASALT << " sophus@" << ORACLE_SOPHUS << " eigen@" << EIGEN_WORLD_VERSION
       << '.' << EIGEN_MAJOR_VERSION << '.' << EIGEN_MINOR_VERSION << " gcc@" << __GNUC__ << '.' << __GNUC_MINOR__
       << '.' << __GNUC_PATCHLEVEL__ << " glibc@" << gnu_get_libc_version();
  return text.str();
}

int run(const std::vector<std::string>& args) {
  if (args.size() == 1 && args[0] == "--version") {
    std::cout << pin() << '\n';
    return std::cout.flush() ? 0 : 1;
  }
  if (args.size() < 3 || args[0] != "--out") {
    std::cerr << "usage: basalt_runner --version | --out DIR FILE.jsonl…\n";
    return 2;
  }
  const std::string dir = args[1];
  std::size_t answered = 0, skipped = 0;
  Counts counts;
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
    answer_file(*id, path, dir + "/" + stem + ".jsonl", counts);
    ++answered;
  }
  std::cerr << "answered " << answered << " ids in binary64; " << counts.off_grid
            << " records off the integer-nanosecond grid and " << counts.refused << " refused (nan); " << skipped
            << " files of other ids skipped\n";
  return 0;
}

}  // namespace
}  // namespace oracle

int main(int argc, char** argv) {
  try {
    return oracle::run(std::vector<std::string>(argv + 1, argv + argc));
  } catch (const std::exception& e) {
    std::cerr << "basalt_runner: " << e.what() << '\n';
    return 1;
  }
}
