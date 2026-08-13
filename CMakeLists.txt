cmake_minimum_required(VERSION 3.28)
# VERSION feeds cramConfigVersion.cmake, so bumping it is what lets a consumer
# write find_package(cram 0.1) and have the request actually checked. See the
# compatibility note where that file is generated.
project(cram_depletion VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Defines CMAKE_INSTALL_INCLUDEDIR / LIBDIR / BINDIR. Included this early
# because target_include_directories(cram) below needs INCLUDEDIR for its
# INSTALL_INTERFACE, not just the install rules at the end of this file.
include(GNUInstallDirs)

# Cache FetchContent downloads outside the build tree so they survive
# clean builds and can be shared across build directories.
set(FETCHCONTENT_BASE_DIR ${CMAKE_SOURCE_DIR}/external CACHE PATH
    "Directory where FetchContent stores downloaded sources")

# All project options are CRAM_-prefixed so they don't collide with the
# options of packages pulled in via FetchContent (Eigen, googletest, ENDFtk),
# which share CMake's global option/cache namespace.
option(CRAM_WITH_ENDFTK      "Build the ENDFtk-based data reader"        OFF)
option(CRAM_ENABLE_TESTS     "Build the unit tests"                      ON)
option(CRAM_ENABLE_BENCHMARKS "Build the Google Benchmark perf suite"    OFF)
option(CRAM_ENABLE_COVERAGE  "Instrument for code coverage (gcov/gcovr)" OFF)
option(CRAM_ENABLE_SANITIZERS "Build with Address + UB sanitizers"       OFF)
# ThreadSanitizer needs its own build: it cannot be combined with ASan, so it
# gets a separate option (and a separate CI job) rather than joining the one
# above.
option(CRAM_ENABLE_TSAN      "Build with ThreadSanitizer"                OFF)
option(CRAM_ENABLE_CLANG_TIDY "Run clang-tidy during the build"          OFF)
# Off by default so a future compiler inventing a new warning never breaks a
# user's build. CI turns it on, which is where new warnings must be caught.
option(CRAM_WERROR           "Treat compiler warnings as errors"         OFF)
# CRAM_ENABLE_INSTALL is declared further down: its default depends on whether
# Eigen ends up found or fetched, which is not known until after the Eigen
# block below.

set(CRAM_EIGEN_MIN_VERSION "3.4" CACHE STRING
    "Minimum acceptable version of a system-installed Eigen")
set(CRAM_EIGEN_GIT_REPOSITORY "https://gitlab.com/libeigen/eigen.git" CACHE STRING
    "Git repository used when Eigen must be fetched instead of found")
set(CRAM_EIGEN_GIT_TAG "5.0.1" CACHE STRING
    "Git tag, branch, or commit of Eigen used when it must be fetched")

# --- Eigen ------------------------------------------------------------------
# Eigen's package version file is generated with SameMajorVersion
# compatibility, so a find_package(Eigen3 3.4) request rejects a system
# Eigen 5.x outright. Search unconstrained and apply our own floor afterwards.
find_package(Eigen3 QUIET NO_MODULE)
if(Eigen3_FOUND AND Eigen3_VERSION VERSION_LESS CRAM_EIGEN_MIN_VERSION)
  message(STATUS "Ignoring system Eigen ${Eigen3_VERSION} "
                 "(< ${CRAM_EIGEN_MIN_VERSION})")
  set(Eigen3_FOUND FALSE)
endif()

if(Eigen3_FOUND)
  message(STATUS "Using system Eigen ${Eigen3_VERSION}")
  set(CRAM_EIGEN_FETCHED FALSE)
else()
  # Recorded because it decides whether cram can be installed at all: a fetched
  # Eigen is a build-tree target, not an installed package, so an exported cram
  # naming it is rejected by install(EXPORT). See the install policy block
  # below, which reads this to pick CRAM_ENABLE_INSTALL's default.
  set(CRAM_EIGEN_FETCHED TRUE)
  message(STATUS "Fetching Eigen ${CRAM_EIGEN_GIT_TAG} from "
                 "${CRAM_EIGEN_GIT_REPOSITORY}")
  include(FetchContent)
  # Eigen 3.4 declares an unconditional option(BUILD_TESTING ... ON) and then
  # registers its ~900-test suite into our CTest project as targets we never
  # build, so `ctest` would try to run nonexistent binaries. Eigen >= 5 gates
  # that (and BLAS/LAPACK/docs/demos) on PROJECT_IS_TOP_LEVEL and needs no
  # help from us, but 3.4-era tags are still reachable via CRAM_EIGEN_GIT_TAG,
  # so keep forcing BUILD_TESTING off in the cache for Eigen's configure (a
  # cache entry — which option() never overrides — is required because Eigen
  # 3.4's old cmake_minimum_required resets CMP0077). Restore the user's own
  # value afterwards so include(CTest) below still sees it.
  if(DEFINED CACHE{BUILD_TESTING})
    set(_cram_saved_build_testing "$CACHE{BUILD_TESTING}")
  endif()
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

  # Eigen 3.4 calls export(PACKAGE Eigen3) unconditionally, and its old
  # cmake_minimum_required leaves CMP0090 at OLD, so the fetched build tree
  # gets written into the user-wide CMake package registry. Every later
  # find_package(Eigen3) on that machine — in this project or any other —
  # would then resolve to a throwaway build directory. Eigen >= 5 guards the
  # export on PROJECT_IS_TOP_LEVEL; suppress it for the tags that don't.
  set(CMAKE_EXPORT_NO_PACKAGE_REGISTRY ON)

  # GIT_SHALLOW can only fetch a named branch or tag, not a bare commit.
  if(CRAM_EIGEN_GIT_TAG MATCHES "^[0-9a-fA-F]{7,40}$")
    set(_cram_eigen_shallow FALSE)
  else()
    set(_cram_eigen_shallow TRUE)
  endif()

  FetchContent_Declare(eigen
    GIT_REPOSITORY ${CRAM_EIGEN_GIT_REPOSITORY}
    GIT_TAG        ${CRAM_EIGEN_GIT_TAG}
    GIT_SHALLOW    ${_cram_eigen_shallow}
    EXCLUDE_FROM_ALL
    SYSTEM)
  FetchContent_MakeAvailable(eigen)

  if(DEFINED _cram_saved_build_testing)
    set(BUILD_TESTING "${_cram_saved_build_testing}" CACHE BOOL
        "Build the testing tree." FORCE)
  else()
    unset(BUILD_TESTING CACHE)
  endif()
endif()

# --- install policy ---------------------------------------------------------
# Two independent conditions have to hold before install rules make sense.
#
# Top-level: when cram is pulled into a larger project via add_subdirectory or
# FetchContent, its install rules would inject cram's headers and package config
# into that project's install tree, which is almost never wanted.
#
# No fetched dependencies: a fetched dependency is a target in *this* build tree
# and is never installed, so exporting cram -- whose PUBLIC interface names
# Eigen3::Eigen (really the `eigen` target) and optionally njoy::ENDFtk -- is not
# merely non-relocatable but a hard install(EXPORT) error:
#
#   install(EXPORT "cramTargets" ...) includes target "cram" which requires
#   target "eigen" that is not in any export set.
#
# Defaulting off in that case keeps `cmake -S . -B build` working for the common
# developer setup with no system Eigen, instead of failing the configure over an
# install nobody asked for. Forcing it on anyway is diagnosed in the install
# section rather than left to surface as the error above.
if(PROJECT_IS_TOP_LEVEL AND NOT CRAM_EIGEN_FETCHED AND NOT CRAM_WITH_ENDFTK)
  set(_cram_install_default ON)
else()
  set(_cram_install_default OFF)
endif()
option(CRAM_ENABLE_INSTALL "Generate install and export rules"
       ${_cram_install_default})

# --- shared instrumentation flags (coverage / sanitizers) ------------------
add_library(cram_quality_flags INTERFACE)
if(CRAM_ENABLE_COVERAGE)
  target_compile_options(cram_quality_flags INTERFACE --coverage -O0 -g)
  target_link_options(cram_quality_flags INTERFACE --coverage)
endif()
if(CRAM_ENABLE_SANITIZERS)
  target_compile_options(cram_quality_flags INTERFACE
    -fsanitize=address,undefined -fno-omit-frame-pointer -g)
  target_link_options(cram_quality_flags INTERFACE -fsanitize=address,undefined)
endif()
if(CRAM_ENABLE_TSAN)
  if(CRAM_ENABLE_SANITIZERS)
    message(FATAL_ERROR
            "CRAM_ENABLE_TSAN and CRAM_ENABLE_SANITIZERS are mutually exclusive: "
            "ThreadSanitizer cannot be combined with AddressSanitizer.")
  endif()
  target_compile_options(cram_quality_flags INTERFACE
    -fsanitize=thread -fno-omit-frame-pointer -g)
  target_link_options(cram_quality_flags INTERFACE -fsanitize=thread)
endif()

# --- warning flags ----------------------------------------------------------
# Linked PRIVATE into every first-party target -- library, tests, driver and
# benchmarks alike. Previously these were set directly on `cram` and were
# therefore PRIVATE to it, so nothing in tests/, cram-apps/ or bench/ was ever
# compiled with -Wall: a warning in test code could not be caught by any CI job.
# Keep the linkage PRIVATE so consumers of `cram` do not inherit our warning
# policy.
add_library(cram_warning_flags INTERFACE)
if(MSVC)
  target_compile_options(cram_warning_flags INTERFACE /W4)
  if(CRAM_WERROR)
    target_compile_options(cram_warning_flags INTERFACE /WX)
  endif()
else()
  target_compile_options(cram_warning_flags INTERFACE -Wall -Wextra)
  if(CRAM_WERROR)
    target_compile_options(cram_warning_flags INTERFACE -Werror)
  endif()
endif()

# --- core library -----------------------------------------------------------
add_library(cram
  cram/cram.cpp
  cram/cram_solver.cpp
  cram/burnup_matrix.cpp
  cram/endf_reader.cpp)
# Headers are included as "cram/<name>.hpp", so the include root is the
# directory *containing* cram/ -- the source root in the build tree, and the
# install includedir once installed. A bare ${CMAKE_CURRENT_SOURCE_DIR} here
# is not merely wrong for consumers but fatal: install(EXPORT) rejects an
# absolute path into the source tree outright.
target_include_directories(cram PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(cram PUBLIC Eigen3::Eigen cram_quality_flags)
# BUILD_INTERFACE because our warning policy is a build-time concern that has no
# business in the installed package. PRIVATE alone is not enough to keep it out:
# cram is a static library, so CMake records its PRIVATE dependencies in
# INTERFACE_LINK_LIBRARIES as $<LINK_ONLY:cram_warning_flags>, and install(EXPORT)
# then rejects the export because that target is not in the set. Adding it to the
# set would satisfy CMake but ship consumers a cram::cram_warning_flags carrying
# nothing usable -- LINK_ONLY drops the compile options that are its whole point.
target_link_libraries(cram PRIVATE $<BUILD_INTERFACE:cram_warning_flags>)
target_compile_features(cram PUBLIC cxx_std_20)

if(CRAM_ENABLE_CLANG_TIDY)
  find_program(CLANG_TIDY_EXE NAMES clang-tidy)
  if(CLANG_TIDY_EXE)
    set_target_properties(cram PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
  endif()
endif()

if(CRAM_WITH_ENDFTK)
  include(FetchContent)
  # Force spdlog (pulled in transitively by njoy/tools) to use std::format
  # instead of its bundled fmt v10, which trips constexpr format-string
  # checks of recent compilers. Set as cache before MakeAvailable so spdlog
  # itself is built with it — propagates to consumers via spdlog::spdlog.
  set(SPDLOG_USE_STD_FORMAT ON CACHE BOOL "" FORCE)
  # SYSTEM, as for Eigen, googletest and benchmark. endf_reader.cpp includes
  # <ENDFtk.hpp>, so without it every diagnostic raised inside ENDFtk's headers
  # (and its njoy/tools dependency) is attributed to our translation unit: with
  # -Werror on, ENDFtk v1.2.0 alone breaks the build with sign-compare, reorder,
  # pessimizing-move and unused-parameter errors from code we do not own.
  FetchContent_Declare(ENDFtk
    GIT_REPOSITORY https://github.com/njoy/ENDFtk.git
    GIT_TAG        v1.2.0
    GIT_SHALLOW    TRUE
    SYSTEM)
  FetchContent_MakeAvailable(ENDFtk)
  target_compile_definitions(cram PUBLIC WITH_ENDFTK)
  target_link_libraries(cram PUBLIC njoy::ENDFtk)
endif()

# --- driver -----------------------------------------------------------------
add_executable(deplete cram-apps/deplete.cpp)
target_link_libraries(deplete PRIVATE cram cram_warning_flags)

# --- tests ------------------------------------------------------------------
include(CTest)
if(CRAM_ENABLE_TESTS)
  find_package(GTest QUIET)
  if(NOT GTest_FOUND)
    include(FetchContent)
    # SYSTEM matches how Eigen and benchmark are fetched. Without it, gtest's
    # headers are treated as ordinary includes and diagnostics raised inside its
    # macro expansions are attributed to our test files -- which is why a warning
    # could appear on a developer machine that fetches gtest and never appear in
    # CI, where libgtest-dev is installed and resolves under /usr/include.
    FetchContent_Declare(googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG        v1.15.2
      SYSTEM)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
  endif()
  add_subdirectory(tests)
endif()

# --- benchmarks -------------------------------------------------------------
if(CRAM_ENABLE_BENCHMARKS)
  add_subdirectory(bench)
endif()

# --- coverage report target -------------------------------------------------
# Uses gcov + cmake/gcov_to_lcov.py to produce an LCOV tracefile that
# VS Code's CMake Tools extension (cmake.coverageInfoFiles) consumes.
if(CRAM_ENABLE_COVERAGE)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  find_program(GCOV_EXE NAMES gcov REQUIRED)

  set(CRAM_COVERAGE_INFO_FILE "${CMAKE_BINARY_DIR}/coverage/coverage.info")

  # Include test sources as gcov inputs so coverage of header-only inline
  # functions (which only get instantiated in test TUs) is captured. The
  # tests/ and cram-apps/ paths are then --exclude'd from the report so
  # only library code (cram/) shows up in coverage numbers.
  get_target_property(_cram_sources cram SOURCES)
  file(GLOB_RECURSE _test_sources CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")
  set(_coverage_source_args "")
  foreach(_src IN LISTS _cram_sources _test_sources)
    if(NOT IS_ABSOLUTE "${_src}")
      set(_src "${CMAKE_CURRENT_SOURCE_DIR}/${_src}")
    endif()
    list(APPEND _coverage_source_args "--source" "${_src}")
  endforeach()

  add_custom_target(coverage-reset
    COMMAND "${CMAKE_COMMAND}"
      -DCRAM_BUILD_DIR=${CMAKE_BINARY_DIR}
      -DCRAM_COVERAGE_INFO_FILE=${CRAM_COVERAGE_INFO_FILE}
      -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ResetCoverage.cmake"
    COMMENT "Removing stale coverage artifacts"
    VERBATIM)

  add_custom_target(coverage
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${CMAKE_BINARY_DIR}/coverage"
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/gcov_to_lcov.py"
      --build-dir "${CMAKE_BINARY_DIR}"
      --source-dir "${CMAKE_CURRENT_SOURCE_DIR}"
      --object-dir "${CMAKE_BINARY_DIR}"
      --output "${CRAM_COVERAGE_INFO_FILE}"
      "--gcov=${GCOV_EXE}"
      --exclude "${FETCHCONTENT_BASE_DIR}"
      --exclude "${CMAKE_CURRENT_SOURCE_DIR}/tests"
      --exclude "${CMAKE_CURRENT_SOURCE_DIR}/cram-apps"
      ${_coverage_source_args}
    BYPRODUCTS "${CRAM_COVERAGE_INFO_FILE}"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating LCOV coverage at ${CRAM_COVERAGE_INFO_FILE}"
    COMMAND_EXPAND_LISTS
    DEPENDS coverage-reset
    VERBATIM)
endif()

# --- install / export -------------------------------------------------------
if(CRAM_ENABLE_INSTALL)
  include(CMakePackageConfigHelpers)

  # Reached only when the user set CRAM_ENABLE_INSTALL=ON explicitly, since the
  # default above is already off in both cases. Diagnosed here with the remedy
  # rather than left to install(EXPORT)'s "target ... not in any export set",
  # which names the dependency's internal target and not the cause.
  if(CRAM_EIGEN_FETCHED)
    message(FATAL_ERROR
      "CRAM_ENABLE_INSTALL is ON but Eigen was fetched rather than found. A "
      "fetched Eigen is a target in this build tree and is never installed, so "
      "exporting cram against it cannot produce a usable package. Install a "
      "system Eigen >= ${CRAM_EIGEN_MIN_VERSION} and reconfigure, or build with "
      "-DCRAM_ENABLE_INSTALL=OFF.")
  endif()
  if(CRAM_WITH_ENDFTK)
    message(FATAL_ERROR
      "CRAM_ENABLE_INSTALL and CRAM_WITH_ENDFTK are both ON. ENDFtk is always "
      "fetched and never installed, so njoy::ENDFtk cannot appear in an "
      "exported cram. Installing an ENDFtk-enabled cram requires ENDFtk to be "
      "packaged and found via find_package first.")
  endif()

  # cram_quality_flags is a PUBLIC dependency of cram, so it is part of cram's
  # exported interface and has to travel with it -- an export set that names a
  # target it does not contain is a configure-time error. It carries nothing at
  # all unless CRAM_ENABLE_COVERAGE or CRAM_ENABLE_SANITIZERS is on, and in
  # those builds consumers genuinely do need the flags: linking a --coverage
  # static library without --coverage leaves the gcov symbols undefined.
  install(TARGETS cram cram_quality_flags
    EXPORT cramTargets
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

  # The driver is installed but deliberately kept out of the export set: it is
  # an application, and exporting it would put a cram::deplete target into every
  # consumer's namespace for no reason.
  install(TARGETS deplete RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  # Listed rather than globbed so that adding a header is a deliberate decision
  # about whether it is public. cram_poles_internal.hpp is intentionally absent:
  # it is included only by cram/*.cpp, never by an installed header.
  install(FILES
    cram/cram.hpp
    cram/cram_solver.hpp
    cram/chain.hpp
    cram/nuclide.hpp
    cram/endf_reader.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cram)

  install(EXPORT cramTargets
    FILE cramTargets.cmake
    NAMESPACE cram::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cram)

  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cramConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cramConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cram)

  # SameMinorVersion, not the more usual SameMajorVersion: before 1.0 the major
  # is pinned at 0 and it is the minor that moves on a breaking change, so
  # SameMajorVersion would declare every 0.x mutually compatible. This should
  # become SameMajorVersion at 1.0.
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/cramConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion)

  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cramConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/cramConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cram)
endif()

# --- cppcheck static analysis target ---------------------------------------
find_program(CPPCHECK_EXE NAMES cppcheck)
if(CPPCHECK_EXE)
  add_custom_target(cppcheck
    COMMAND ${CPPCHECK_EXE}
            --enable=warning,style,performance,portability
            --std=c++20 --language=c++ --inline-suppr --error-exitcode=1
            -I ${CMAKE_SOURCE_DIR}
            --suppress=missingIncludeSystem
            --suppress=normalCheckLevelMaxBranches
            ${CMAKE_SOURCE_DIR}/cram
    COMMENT "Running cppcheck static analysis")
endif()
