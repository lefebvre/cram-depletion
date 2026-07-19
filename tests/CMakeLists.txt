add_executable(unit_tests
  test_nuclide.cpp
  test_chain.cpp
  test_cram.cpp
  test_cram_solver.cpp)

target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(unit_tests PRIVATE cram cram_quality_flags GTest::gtest_main)

include(GoogleTest)
# PRE_TEST defers test discovery to `ctest` invocation time rather than
# running the freshly built binary as a POST_BUILD step. On Windows
# multi-config generators the POST_BUILD run can misresolve the target
# path and record garbage test entries; PRE_TEST discovers on the same
# machine in the correct configuration.
gtest_discover_tests(unit_tests DISCOVERY_MODE PRE_TEST)

# Integration test for the real ENDFtk reader. Built only with -DWITH_ENDFTK.
if(WITH_ENDFTK)
  add_executable(integration_test integration/integration_test.cpp)
  target_link_libraries(integration_test PRIVATE cram GTest::gtest_main)
  gtest_discover_tests(integration_test
    DISCOVERY_MODE PRE_TEST
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/integration)
endif()
