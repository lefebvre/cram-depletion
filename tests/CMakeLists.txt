add_executable(unit_tests
  test_nuclide.cpp
  test_chain.cpp
  test_cram.cpp
  test_cram_solver.cpp
  test_integrator.cpp
  test_deplete.cpp
  test_chain_xml.cpp)

target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(unit_tests PRIVATE cram cram_quality_flags GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(unit_tests)

# Validation replay against OpenMC-generated VERA pin reference data. The data
# lives under tests/validation/data/ and is produced offline (see
# validation/openmc/); the test SKIPS when it is absent, so this always builds
# and runs without the multi-GB nuclear-data library.
add_executable(validation_test validation/vera_pin_test.cpp)
target_link_libraries(validation_test PRIVATE cram cram_quality_flags GTest::gtest_main)
target_compile_definitions(validation_test PRIVATE
  CRAM_VALIDATION_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/validation/data")
gtest_discover_tests(validation_test)

# Integration test for the real ENDFtk reader. Built only with -DWITH_ENDFTK.
if(WITH_ENDFTK)
  add_executable(integration_test integration/integration_test.cpp)
  target_link_libraries(integration_test PRIVATE cram GTest::gtest_main)
  gtest_discover_tests(integration_test
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/integration)
endif()
