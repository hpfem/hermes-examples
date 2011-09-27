project(capillary-barrier-adapt)
add_executable(${PROJECT_NAME} main.cpp definitions.cpp constitutive_relations.cpp definitions.h)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")