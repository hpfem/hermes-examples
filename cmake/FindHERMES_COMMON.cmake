INCLUDE(FindPackageHandleStandardArgs)

if(WIN64)
  FIND_LIBRARY(HERMES_COMMON_LIBRARY NAMES hermes_common_64 hermes_common_64-debug PATHS ${HERMES_DIRECTORY} /usr/lib /usr/local/lib NO_DEFAULT_PATH)
else()
  FIND_LIBRARY(HERMES_COMMON_LIBRARY NAMES hermes_common hermes_common-debug PATHS ${HERMES_DIRECTORY} /usr/lib /usr/local/lib NO_DEFAULT_PATH)
endif()
  
FIND_PATH(HERMES_COMMON_INCLUDE hermes_common.h ${HERMES_COMMON_INCLUDE_PATH} /usr/include/hermes_common /usr/local/include/hermes_common)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HERMES_COMMON DEFAULT_MSG HERMES_COMMON_LIBRARY HERMES_COMMON_INCLUDE)
