find_package(benchmark QUIET)
if(NOT benchmark_FOUND)
  include(FetchContent)
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG        v1.9.0
    GIT_SHALLOW    TRUE
    EXCLUDE_FROM_ALL
    SYSTEM)
  FetchContent_MakeAvailable(benchmark)
endif()

add_executable(cram_bench bench_cram.cpp)
# testing/ holds the shared cyclic-chain fixture. The benchmarks and the golden
# regression test deliberately measure the *same* matrix structure, so the
# generator lives in one place rather than being duplicated here.
target_include_directories(cram_bench PRIVATE ${CMAKE_SOURCE_DIR}/testing)
target_link_libraries(cram_bench PRIVATE cram cram_warning_flags benchmark::benchmark_main)
