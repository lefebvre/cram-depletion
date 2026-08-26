add_executable(unit_tests
  test_nuclide.cpp
  test_chain.cpp
  test_reaction.cpp
  test_cram.cpp
  test_cram_solver.cpp
  test_cram_golden.cpp
  test_integrator.cpp
  test_deplete.cpp)

# The reader's tests need the reader; the stub-throws test lives in
# test_chain.cpp beside the ENDFtk one and always runs.
if(CRAM_WITH_CHAIN_XML)
  target_sources(unit_tests PRIVATE test_chain_xml.cpp)
endif()

target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                              ${CMAKE_SOURCE_DIR}/testing)
target_link_libraries(unit_tests PRIVATE cram cram_quality_flags cram_warning_flags
                                         GTest::gtest_main)

include(GoogleTest)
# TEST_PREFIX namespaces our tests as `cram.*` in CTest so they're clearly
# distinguished from any test a dependency might register. (Eigen's own suite
# is kept out entirely by forcing BUILD_TESTING off around its configure in
# the top-level CMakeLists.txt.)
# DISCOVERY_MODE PRE_TEST runs the enumeration at test time rather than as a
# POST_BUILD step. Under a sanitizer the difference matters: ThreadSanitizer can
# abort at startup on kernels that hand out more ASLR entropy than its shadow
# mapping allows, and with POST_BUILD discovery that failure breaks the *build*
# rather than the test run — intermittently, since it depends on where the
# loader happens to land. Deferring discovery keeps a sanitizer build from
# failing before anything has been run.
gtest_discover_tests(unit_tests TEST_PREFIX cram. DISCOVERY_MODE PRE_TEST)

# Replays an OpenMC-generated VERA pin depletion (tests/validation/data) and
# checks the engine reproduces OpenMC's number densities. Needs the chain
# reader; SKIPs at run time when the data directory is absent. A separate
# binary so the ~1 MB of reference data is parsed only by the test that uses it.
if(CRAM_WITH_CHAIN_XML)
  add_executable(validation_test validation/vera_pin_test.cpp)
  target_link_libraries(validation_test PRIVATE cram cram_quality_flags cram_warning_flags
                                                GTest::gtest_main)
  target_compile_definitions(validation_test PRIVATE
    CRAM_VALIDATION_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/validation/data")
  gtest_discover_tests(validation_test TEST_PREFIX cram. DISCOVERY_MODE PRE_TEST)
endif()

# Integration test for the real ENDFtk reader. Built only with -DCRAM_WITH_ENDFTK.
if(CRAM_WITH_ENDFTK)
  add_executable(integration_test integration/integration_test.cpp)
  target_link_libraries(integration_test PRIVATE cram cram_warning_flags GTest::gtest_main)
  gtest_discover_tests(integration_test
    TEST_PREFIX cram.
    DISCOVERY_MODE PRE_TEST
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/integration)
endif()
