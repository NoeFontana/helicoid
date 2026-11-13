// The corpus's hex floats (`float.hex()`) in, C's `%a` out: both exact, so no rounding enters a
// record or an answer. `nan`, `inf` and `-inf` carry a non-finite answer to the harness.
#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace oracle {

// `strtod` reads a hex float exactly; anything but a whole `[-]0x…p…` string is an error.
inline double parse_hex(const std::string& s) {
  const bool hex = s.rfind("0x", 0) == 0 || s.rfind("-0x", 0) == 0;
  const bool binary = s.find('p') != std::string::npos;
  char* end = nullptr;
  const double x = std::strtod(s.c_str(), &end);
  if (!hex || !binary || end != s.c_str() + s.size()) {
    throw std::invalid_argument("bad hex float `" + s + "`");
  }
  return x;
}

inline std::string format_hex(double x) {
  if (std::isnan(x)) return "nan";
  if (std::isinf(x)) return x < 0 ? "-inf" : "inf";
  char buffer[64];
  std::snprintf(buffer, sizeof buffer, "%a", x);
  return buffer;
}

}  // namespace oracle
