# Fetch ENDF/B-VIII.0 decay tapes used by the Pusa-2016 actinide
# validation test (Pu-241 -> Am-241 -> Np-237 -> Pa-233 -> U-233).
#
# Hash-pinned downloads from the IAEA Nuclear Data Services public
# distribution. Cached under FETCHCONTENT_BASE_DIR so clean build
# directories don't re-download.

set(PUSA_ACTINIDE_DATA_DIR
    "${FETCHCONTENT_BASE_DIR}/endf_pusa_actinide"
    CACHE INTERNAL "Directory holding the fetched Pusa actinide ENDF tapes")
file(MAKE_DIRECTORY "${PUSA_ACTINIDE_DATA_DIR}")

# (basename, sha256) for IAEA NDS ENDF/B-VIII.0 decay zips.
# URL: https://www-nds.iaea.org/public/download-endf/ENDF-B-VIII.0/decay/<basename>.zip
set(_pusa_tapes
  "decay_3489_91-Pa-233|61608c74ffa64b35b168d49f83b89e5a9ca854bc9e9560f17ecb366c4c2c2fae"
  "decay_3513_92-U-233|3acd582d91e8e0240631463001f84ed66a86400ae637649a252515ee924ffb37"
  "decay_3537_93-Np-237|447cb954fb3336c04b36b867a6d096dd4f983cad9c592748d198e7886630c4fb"
  "decay_3561_94-Pu-241|6a5771d0aed391640a7f9434998dcf106027bcaf847841d85605d79291292b46"
  "decay_3578_95-Am-241|e3739aa638e582d750a582429009e56ccc5e2c476bdbb5d23264cc396cfa4a79"
)

foreach(_entry IN LISTS _pusa_tapes)
  string(REPLACE "|" ";" _parts "${_entry}")
  list(GET _parts 0 _base)
  list(GET _parts 1 _hash)
  set(_zip "${PUSA_ACTINIDE_DATA_DIR}/${_base}.zip")
  set(_dat "${PUSA_ACTINIDE_DATA_DIR}/${_base}.dat")
  if(NOT EXISTS "${_dat}")
    message(STATUS "Fetching ENDF/B-VIII.0 decay tape: ${_base}")
    file(DOWNLOAD
      "https://www-nds.iaea.org/public/download-endf/ENDF-B-VIII.0/decay/${_base}.zip"
      "${_zip}"
      EXPECTED_HASH SHA256=${_hash}
      STATUS _status)
    list(GET _status 0 _err)
    if(NOT _err EQUAL 0)
      list(GET _status 1 _msg)
      message(FATAL_ERROR
        "Failed to fetch ${_base}.zip from IAEA NDS (${_msg}).\n"
        "WITH_ENDFTK=ON requires network access at configure time to download "
        "the Pusa-2016 validation tapes. Re-configure when the network is "
        "reachable, or set -DWITH_ENDFTK=OFF to skip this validation.")
    endif()
    file(ARCHIVE_EXTRACT INPUT "${_zip}" DESTINATION "${PUSA_ACTINIDE_DATA_DIR}")
  endif()
endforeach()
