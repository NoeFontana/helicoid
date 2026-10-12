# syntax=docker/dockerfile:1
# The GTSAM C++ oracle (docs/PHASE1.md §7, docs/decisions/0074). Everything it computes with is
# pinned: the base image by digest, every package by the Ubuntu archive snapshot, GTSAM by commit and
# checksum, compiled here under this file's flags. `cargo xtask conformance --oracle gtsam` builds
# this and runs it offline.
FROM ubuntu@sha256:534baea6a22c03a63003dbc8dbe78fe34bc0d7e595d9a9dc9834884ff530eb55

# The archive as it was at this instant: g++ 13.3.0, cmake 3.28.3, nlohmann-json3-dev 3.11.3,
# libc6 2.39. `ca-certificates` comes first, from the live archive, only to reach the snapshot over
# https.
ARG SNAPSHOT=20261001T000000Z
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates \
 && apt-get update --snapshot "${SNAPSHOT}" \
 && apt-get install -y --no-install-recommends --snapshot "${SNAPSHOT}" \
      g++ cmake make nlohmann-json3-dev \
 && rm -rf /var/lib/apt/lists/*

# GTSAM 4.3.0, its latest release, with the Eigen it bundles (3.4.0). Its options are the defaults
# wherever they reach the numerics (`Rot3` a matrix, `Pose3`'s chart `Expmap`); the rest only drop
# what the oracle never calls (Boost, TBB, METIS, tests, examples, `gtsam_unstable`) and its own
# `-Werror`, which gcc 13's `-Wmaybe-uninitialized` trips in `Gal3.cpp` without Eigen SIMD. No
# contraction into fused multiply-adds and no Eigen SIMD, so the library and the runner are the
# same scalar program on every host.
ARG GTSAM_VERSION=4.3.0
ARG GTSAM_COMMIT=71a25ca36c084cbad1f872e812d6d97fbadfdb05
ARG CXXFLAGS_ORACLE="-ffp-contract=off -DEIGEN_DONT_VECTORIZE"
ADD --checksum=sha256:fa79b78a2cdd2561d4f93b75a61294a39c4d3a2d2e6e70dc51f91e6253966f96 \
    https://github.com/borglab/gtsam/archive/${GTSAM_COMMIT}.tar.gz /tmp/gtsam.tar.gz
RUN set -eu; \
    mkdir -p /tmp/gtsam; \
    tar -xzf /tmp/gtsam.tar.gz -C /tmp/gtsam --strip-components=1; \
    cmake -S /tmp/gtsam -B /tmp/gtsam-build \
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/gtsam \
      -DCMAKE_CXX_FLAGS="${CXXFLAGS_ORACLE}" \
      -DBUILD_SHARED_LIBS=OFF -DGTSAM_BUILD_TYPE_POSTFIXES=OFF \
      -DGTSAM_BUILD_WITH_MARCH_NATIVE=OFF -DGTSAM_BUILD_WITH_CCACHE=OFF -DGTSAM_BUILD_WITH_WERROR=OFF \
      -DGTSAM_USE_BOOST_FEATURES=OFF -DGTSAM_ENABLE_BOOST_SERIALIZATION=OFF \
      -DGTSAM_WITH_TBB=OFF -DGTSAM_SUPPORT_NESTED_DISSECTION=OFF \
      -DGTSAM_BUILD_UNSTABLE=OFF -DGTSAM_BUILD_TESTS=OFF \
      -DGTSAM_BUILD_EXAMPLES_ALWAYS=OFF -DGTSAM_BUILD_TIMING_ALWAYS=OFF \
      -DGTSAM_BUILD_PYTHON=OFF -DGTSAM_INSTALL_MATLAB_TOOLBOX=OFF; \
    cmake --build /tmp/gtsam-build --parallel "$(nproc)"; \
    cmake --install /tmp/gtsam-build; \
    rm -rf /tmp/gtsam /tmp/gtsam-build /tmp/gtsam.tar.gz

# The runner, against the installed library and under the same flags, with warnings denied; the
# conversions' hand-computed tests run first and stop the build on a failure.
COPY CMakeLists.txt convert.hpp hexfloat.hpp runner.cpp convert_test.cpp /src/
RUN set -eu; \
    cmake -S /src -B /tmp/runner-build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH=/opt/gtsam -DCMAKE_CXX_FLAGS="${CXXFLAGS_ORACLE}" \
      -DORACLE_GTSAM="${GTSAM_VERSION}+${GTSAM_COMMIT}"; \
    cmake --build /tmp/runner-build --parallel "$(nproc)"; \
    /tmp/runner-build/convert_test; \
    install -m 0755 /tmp/runner-build/gtsam_runner /usr/local/bin/gtsam_runner; \
    rm -rf /tmp/runner-build

ENTRYPOINT ["/usr/local/bin/gtsam_runner"]
