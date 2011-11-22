project(capillary-barrier-adapt)
add_executable(${PROJECT_NAME} main.cpp definitions.cpp extras.cpp definitions.h)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")