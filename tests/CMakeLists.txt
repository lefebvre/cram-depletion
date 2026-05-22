add_executable(unit_tests
  test_nuclide.cpp
  test_chain.cpp
  test_cram.cpp)

target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(unit_tests PRIVATE cram cram_quality_flags GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(unit_tests)

# Integration test for the real ENDFtk reader. Built only with -DWITH_ENDFTK.
if(WITH_ENDFTK)
  add_executable(integration_test integration/integration_test.cpp)
  target_link_libraries(integration_test PRIVATE cram GTest::gtest_main)
  gtest_discover_tests(integration_test
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/integration)

  # Pusa-2016 actinide validation: reads tapes fetched into
  # ${PUSA_ACTINIDE_DATA_DIR} by cmake/FetchEndfActinides.cmake.
  add_executable(pusa_actinide_test integration/test_pusa_2016_actinide.cpp)
  target_link_libraries(pusa_actinide_test PRIVATE cram GTest::gtest_main)
  gtest_discover_tests(pusa_actinide_test
    WORKING_DIRECTORY ${PUSA_ACTINIDE_DATA_DIR})
endif()
