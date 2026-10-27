add_executable(unit_tests
  test_nuclide.cpp
  test_chain.cpp
  test_reaction.cpp
  test_cram.cpp
  test_cram_poles.cpp
  test_cram_solver.cpp
  test_cram_golden.cpp
  test_integrator.cpp
  test_deplete.cpp
  test_adjoint.cpp
  test_sensitivity.cpp)

# The reader's tests need the reader. Its counterpart, ChainXml.StubsThrowWithoutReader,
# lives in test_chain.cpp beside the ENDFtk one under #ifndef CRAM_WITH_CHAIN_XML,
# so exactly one of the two is compiled in any configuration: with the reader,
# the tests that parse; without it, the one that asserts the loaders throw
# rather than leaving the chain silently empty.
if(CRAM_WITH_CHAIN_XML)
  target_sources(unit_tests PRIVATE test_chain_xml.cpp)
endif()

target_include_directories(unit_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                              ${PROJECT_SOURCE_DIR}/testing)
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

# Compile-fail tests for the structs whose required fields have no default
# member initializer: compile_fail/missing_field.cpp explains the guarantee.
# Each case builds one target from inside ctest and passes only if the build
# reports -Wmissing-field-initializers (or Clang's designated-initializer
# subgroup of it) -- matching the diagnostic, not just a failed exit code, so an
# unrelated compile error cannot pass for one. The closing bracket is part of
# the match because the flag alone also appears in the command line a failed
# Ninja build echoes. The control build names every
# field and must succeed, which keeps the snippet itself honest.
#
# GCC and Clang only: MSVC has no warning for an omitted aggregate member, so
# on MSVC the guarantee does not hold at all and there is nothing to test.
if(NOT MSVC)
  set(_cram_missing_field_cases Zai DecayMode DecayData FissionYields ReactionXS ChainReaction)
  add_executable(compile_fail_control EXCLUDE_FROM_ALL compile_fail/missing_field.cpp)
  target_link_libraries(compile_fail_control PRIVATE cram)
  target_compile_options(compile_fail_control PRIVATE -Werror=missing-field-initializers)
  add_test(NAME cram.CompileFail.ControlBuilds
           COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target compile_fail_control
                   --config $<CONFIG>)
  set(_cram_omit 0)
  foreach(_case IN LISTS _cram_missing_field_cases)
    math(EXPR _cram_omit "${_cram_omit} + 1")
    add_executable(compile_fail_${_case} EXCLUDE_FROM_ALL compile_fail/missing_field.cpp)
    target_link_libraries(compile_fail_${_case} PRIVATE cram)
    target_compile_definitions(compile_fail_${_case} PRIVATE CRAM_OMIT=${_cram_omit})
    target_compile_options(compile_fail_${_case} PRIVATE -Werror=missing-field-initializers)
    add_test(NAME cram.CompileFail.MissingFieldIn${_case}
             COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target compile_fail_${_case}
                     --config $<CONFIG>)
    set_tests_properties(cram.CompileFail.MissingFieldIn${_case} PROPERTIES
      PASS_REGULAR_EXPRESSION "missing-(designated-)?field-initializers\]")
  endforeach()
endif()

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
