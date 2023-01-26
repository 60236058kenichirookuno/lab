#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "LBFGS::lib" for configuration "Debug"
set_property(TARGET LBFGS::lib APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(LBFGS::lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/lbfgsd.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS LBFGS::lib )
list(APPEND _IMPORT_CHECK_FILES_FOR_LBFGS::lib "${_IMPORT_PREFIX}/lib/lbfgsd.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
