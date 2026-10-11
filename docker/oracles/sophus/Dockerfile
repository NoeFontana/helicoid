# syntax=docker/dockerfile:1
# The Sophus C++ oracle (docs/PHASE1.md §7, docs/decisions/0072). Everything it computes with is
# pinned: the base image by digest, every package by the Ubuntu archive snapshot, Sophus by commit and
# checksum. `cargo xtask conformance --oracle sophus` builds this and runs it offline.
FROM ubuntu@sha256:534baea6a22c03a63003dbc8dbe78fe34bc0d7e595d9a9dc9834884ff530eb55

# The archive as it was at this instant: g++ 13.3.0, libeigen3-dev 3.4.0, nlohmann-json3-dev
# 3.11.3, libc6 2.39. `ca-certificates` comes first, from the live archive, only to reach the
# snapshot over https.
ARG SNAPSHOT=20261001T000000Z
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates \
 && apt-get update --snapshot "${SNAPSHOT}" \
 && apt-get install -y --no-install-recommends --snapshot "${SNAPSHOT}" \
      g++ libeigen3-dev nlohmann-json3-dev \
 && rm -rf /var/lib/apt/lists/*

# Sophus 1.24.6, its latest release; header-only, so only its headers are kept.
ARG SOPHUS_VERSION=1.24.6
ARG SOPHUS_COMMIT=d0b7315a0d90fc6143defa54596a3a95d9fa10ec
ADD --checksum=sha256:876e24fd402ac4116cd7863685225388b77be86ec356b777bf68defd79decc90 \
    https://github.com/strasdat/Sophus/archive/${SOPHUS_COMMIT}.tar.gz /tmp/sophus.tar.gz
RUN mkdir -p /opt/sophus \
 && tar -xzf /tmp/sophus.tar.gz -C /opt/sophus --strip-components=1 "Sophus-${SOPHUS_COMMIT}/sophus" \
 && rm /tmp/sophus.tar.gz

# No contraction into fused multiply-adds and no Eigen SIMD, so the oracle is the same scalar
# program on every host; a failed `SOPHUS_ENSURE` throws (`refused.cpp`) instead of aborting.
COPY convert.hpp hexfloat.hpp refused.cpp runner.cpp convert_test.cpp /src/
RUN set -eu; \
    flags="-std=c++17 -O2 -Wall -Wextra -Wpedantic -Werror -ffp-contract=off \
      -DEIGEN_DONT_VECTORIZE -DSOPHUS_ENABLE_ENSURE_HANDLER \
      -isystem /usr/include/eigen3 -isystem /opt/sophus"; \
    g++ $flags -o /usr/local/bin/convert_test /src/convert_test.cpp /src/refused.cpp; \
    /usr/local/bin/convert_test; \
    g++ $flags -DORACLE_SOPHUS="\"${SOPHUS_VERSION}+${SOPHUS_COMMIT}\"" \
      -o /usr/local/bin/sophus_runner /src/runner.cpp /src/refused.cpp

ENTRYPOINT ["/usr/local/bin/sophus_runner"]
