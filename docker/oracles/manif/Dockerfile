# syntax=docker/dockerfile:1
# The manif oracle (docs/PHASE1.md §7, docs/decisions/0075). Everything it computes with is pinned:
# the base image by digest, every package by the Ubuntu archive snapshot, manif by commit and
# checksum. `cargo xtask conformance --oracle manif` builds this and runs it offline.
FROM ubuntu@sha256:534baea6a22c03a63003dbc8dbe78fe34bc0d7e595d9a9dc9834884ff530eb55

# The archive as it was at this instant: g++ 13.3.0, libeigen3-dev 3.4.0, nlohmann-json3-dev
# 3.11.3, libc6 2.39, as the Sophus image's (docs/decisions/0072). `ca-certificates` comes first,
# from the live archive, only to reach the snapshot over https.
ARG SNAPSHOT=20261001T000000Z
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates \
 && apt-get update --snapshot "${SNAPSHOT}" \
 && apt-get install -y --no-install-recommends --snapshot "${SNAPSHOT}" \
      g++ libeigen3-dev nlohmann-json3-dev \
 && rm -rf /var/lib/apt/lists/*

# manif 0.0.5, its latest release; header-only, so only its headers and the `tl::optional` it
# vendors are kept.
ARG MANIF_VERSION=0.0.5
ARG MANIF_COMMIT=c49fec8150e62821d9be7e2a3ac49db75bf7bb92
ADD --checksum=sha256:5118fba84f286c4172f0820c9980d732187c471c44238db0243e376e9bd7b6b6 \
    https://github.com/artivis/manif/archive/${MANIF_COMMIT}.tar.gz /tmp/manif.tar.gz
RUN mkdir -p /opt/manif \
 && tar -xzf /tmp/manif.tar.gz -C /opt/manif --strip-components=1 \
      "manif-${MANIF_COMMIT}/include" "manif-${MANIF_COMMIT}/external/tl" \
 && rm /tmp/manif.tar.gz

# No contraction into fused multiply-adds and no Eigen SIMD, so the oracle is the same scalar
# program on every host; no `NDEBUG`, so manif's unit-norm checks throw (`convert.hpp`).
COPY convert.hpp hexfloat.hpp runner.cpp convert_test.cpp /src/
RUN set -eu; \
    flags="-std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -ffp-contract=off -DEIGEN_DONT_VECTORIZE \
      -isystem /usr/include/eigen3 -isystem /opt/manif/include -isystem /opt/manif/external/tl"; \
    g++ $flags -o /usr/local/bin/convert_test /src/convert_test.cpp; \
    /usr/local/bin/convert_test; \
    g++ $flags -DORACLE_MANIF="\"${MANIF_VERSION}+${MANIF_COMMIT}\"" \
      -o /usr/local/bin/manif_runner /src/runner.cpp

ENTRYPOINT ["/usr/local/bin/manif_runner"]
