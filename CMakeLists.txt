cmake_minimum_required(VERSION 3.28)
project(cram_depletion LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Cache FetchContent downloads outside the build tree so they survive
# clean builds and can be shared across build directories.
set(FETCHCONTENT_BASE_DIR ${CMAKE_SOURCE_DIR}/external CACHE PATH
    "Directory where FetchContent stores downloaded sources")

find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  include(FetchContent)
  FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    GIT_SHALLOW    TRUE
    EXCLUDE_FROM_ALL
    SYSTEM)
  FetchContent_MakeAvailable(eigen)
endif()

option(WITH_ENDFTK      "Build the ENDFtk-based data reader"        OFF)
option(ENABLE_TESTS     "Build the unit tests"                      ON)
option(ENABLE_BENCHMARKS "Build the Google Benchmark perf suite"    OFF)
option(ENABLE_COVERAGE  "Instrument for code coverage (gcov/gcovr)" OFF)
option(ENABLE_SANITIZERS "Build with Address + UB sanitizers"       OFF)
option(ENABLE_CLANG_TIDY "Run clang-tidy during the build"          OFF)

# --- shared instrumentation flags (coverage / sanitizers) ------------------
add_library(cram_quality_flags INTERFACE)
if(ENABLE_COVERAGE)
  target_compile_options(cram_quality_flags INTERFACE --coverage -O0 -g)
  target_link_options(cram_quality_flags INTERFACE --coverage)
endif()
if(ENABLE_SANITIZERS)
  target_compile_options(cram_quality_flags INTERFACE
    -fsanitize=address,undefined -fno-omit-frame-pointer -g)
  target_link_options(cram_quality_flags INTERFACE -fsanitize=address,undefined)
endif()

# --- core library -----------------------------------------------------------
add_library(cram
  cram/cram.cpp
  cram/burnup_matrix.cpp
  cram/endf_reader.cpp)
target_include_directories(cram PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(cram PUBLIC Eigen3::Eigen cram_quality_flags)
target_compile_features(cram PUBLIC cxx_std_20)
target_compile_options(cram PRIVATE -Wall -Wextra)

if(ENABLE_CLANG_TIDY)
  find_program(CLANG_TIDY_EXE NAMES clang-tidy)
  if(CLANG_TIDY_EXE)
    set_target_properties(cram PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
  endif()
endif()

if(WITH_ENDFTK)
  include(FetchContent)
  # Force spdlog (pulled in transitively by njoy/tools) to use std::format
  # instead of its bundled fmt v10, which trips constexpr format-string
  # checks of recent compilers. Set as cache before MakeAvailable so spdlog
  # itself is built with it — propagates to consumers via spdlog::spdlog.
  set(SPDLOG_USE_STD_FORMAT ON CACHE BOOL "" FORCE)
  FetchContent_Declare(ENDFtk
    GIT_REPOSITORY https://github.com/njoy/ENDFtk.git
    GIT_TAG        v1.2.0
    GIT_SHALLOW    TRUE)
  FetchContent_MakeAvailable(ENDFtk)
  target_compile_definitions(cram PUBLIC WITH_ENDFTK)
  target_link_libraries(cram PUBLIC njoy::ENDFtk)
endif()

# --- driver -----------------------------------------------------------------
add_executable(deplete cram-apps/deplete.cpp)
target_link_libraries(deplete PRIVATE cram)

# --- tests ------------------------------------------------------------------
include(CTest)
if(ENABLE_TESTS)
  find_package(GTest QUIET)
  if(NOT GTest_FOUND)
    include(FetchContent)
    FetchContent_Declare(googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG        v1.15.2)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
  endif()
  add_subdirectory(tests)
endif()

# --- benchmarks -------------------------------------------------------------
if(ENABLE_BENCHMARKS)
  add_subdirectory(bench)
endif()

# --- coverage report target -------------------------------------------------
# Uses gcov + cmake/gcov_to_lcov.py to produce an LCOV tracefile that
# VS Code's CMake Tools extension (cmake.coverageInfoFiles) consumes.
if(ENABLE_COVERAGE)
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
