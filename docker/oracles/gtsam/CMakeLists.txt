# The GTSAM oracle's runner and its conversion tests (docs/decisions/0074), built against the
# installed static `gtsam` so its public definitions reach this code as they reached the library.
cmake_minimum_required(VERSION 3.16)
project(gtsam_oracle CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(GTSAM REQUIRED CONFIG)
find_package(nlohmann_json REQUIRED CONFIG)

set(ORACLE_GTSAM "" CACHE STRING "the pinned GTSAM version and commit")
if(ORACLE_GTSAM STREQUAL "")
  message(FATAL_ERROR "ORACLE_GTSAM is set by the Dockerfile")
endif()

foreach(target convert_test gtsam_runner)
  add_executable(${target})
  target_link_libraries(${target} PRIVATE gtsam nlohmann_json::nlohmann_json)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
endforeach()
target_sources(convert_test PRIVATE convert_test.cpp)
target_sources(gtsam_runner PRIVATE runner.cpp)
target_compile_definitions(gtsam_runner PRIVATE "ORACLE_GTSAM=\"${ORACLE_GTSAM}\"")
