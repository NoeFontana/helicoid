// `SOPHUS_ENABLE_ENSURE_HANDLER`'s hook: a failed `SOPHUS_ENSURE` is a refusal of the record, which
// the runner answers with `nan` (docs/decisions/0072), instead of Sophus's default `abort`.
#include "convert.hpp"

namespace Sophus {
void ensureFailed(char const* function, char const* file, int line, char const* description) {
  throw oracle::Refused(std::string(function) + " (" + file + ":" + std::to_string(line) +
                        "): " + description);
}
}  // namespace Sophus
