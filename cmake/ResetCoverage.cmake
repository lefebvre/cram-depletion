if(NOT DEFINED CRAM_BUILD_DIR)
  message(FATAL_ERROR "CRAM_BUILD_DIR is required")
endif()
if(NOT DEFINED CRAM_COVERAGE_INFO_FILE)
  message(FATAL_ERROR "CRAM_COVERAGE_INFO_FILE is required")
endif()

file(GLOB_RECURSE _stale
  LIST_DIRECTORIES FALSE
  "${CRAM_BUILD_DIR}/*.gcda"
  "${CRAM_BUILD_DIR}/*.gcov"
  "${CRAM_BUILD_DIR}/*.gcov.json.gz")

if(_stale)
  file(REMOVE ${_stale})
endif()
file(REMOVE "${CRAM_COVERAGE_INFO_FILE}")
