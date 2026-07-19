add_executable(unit_tests
  test_nuclide.cpp
  test_chain.cpp
  test_cram.cpp
  test_cram_solver.cpp)

target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(unit_tests PRIVATE cram cram_quality_flags GTest::gtest_main)

include(GoogleTest)
# TEST_PREFIX namespaces our tests as `cram.*` in CTest so they're clearly
# distinguished from any test a dependency might register. (Eigen's own suite
# is kept out entirely by forcing BUILD_TESTING off around its configure in
# the top-level CMakeLists.txt.)
gtest_discover_tests(unit_tests TEST_PREFIX cram.)

# Integration test for the real ENDFtk reader. Built only with -DCRAM_WITH_ENDFTK.
if(CRAM_WITH_ENDFTK)
  add_executable(integration_test integration/integration_test.cpp)
  target_link_libraries(integration_test PRIVATE cram GTest::gtest_main)
  gtest_discover_tests(integration_test
    TEST_PREFIX cram.
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/integration)
endif()
