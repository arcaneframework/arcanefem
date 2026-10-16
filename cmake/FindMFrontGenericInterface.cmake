# - Try to find the MGIS (MFrontGenericInterfaceSupport) library.
#
# Once done this module will define:
#
#   MFrontGenericInterface_FOUND - whether MGIS was found. The imported
#   target mgis::MFrontGenericInterface is defined on success.
#
# The following variables can be used as hints:
#
#   MFrontGenericInterface_DIR or MGIS_DIR - path of the MGIS installation
#   prefix, or of the directory containing MFrontGenericInterfaceConfig.cmake
#   (usually <prefix>/share/mgis/cmake).

find_path(MFrontGenericInterface_CONFIG_FILE
  NAMES MFrontGenericInterfaceConfig.cmake
  HINTS ${MFrontGenericInterface_DIR} $ENV{MFrontGenericInterface_DIR}
        ${MGIS_DIR} $ENV{MGIS_DIR}
  PATH_SUFFIXES share/mgis/cmake
  DOC "The MFrontGenericInterface configuration file")
mark_as_advanced(MFrontGenericInterface_CONFIG_FILE)

set(MFrontGenericInterface_FOUND FALSE)

if (MFrontGenericInterface_CONFIG_FILE)
  include("${MFrontGenericInterface_CONFIG_FILE}/MFrontGenericInterfaceConfig.cmake")
  if (TARGET mgis::MFrontGenericInterface)
    set(MFrontGenericInterface_FOUND TRUE)
    set(MFrontGenericInterface_DIR "${MFrontGenericInterface_CONFIG_FILE}" CACHE PATH
      "The directory containing the MFrontGenericInterface package configuration" FORCE)
  endif ()
endif ()

if (NOT MFrontGenericInterface_FOUND)
  set(MFrontGenericInterface_FAIL_MESSAGE
    "MFrontGenericInterface (MGIS) could not be found. Please set MGIS_DIR to the MGIS installation prefix")
  if (MFrontGenericInterface_FIND_REQUIRED)
    message(FATAL_ERROR ${MFrontGenericInterface_FAIL_MESSAGE})
  elseif (NOT MFrontGenericInterface_FIND_QUIETLY)
    message(STATUS ${MFrontGenericInterface_FAIL_MESSAGE})
  endif ()
endif ()
