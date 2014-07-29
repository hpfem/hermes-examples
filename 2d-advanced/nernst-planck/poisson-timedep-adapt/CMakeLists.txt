project(nernst-planck-poisson-timedep-adapt)

add_executable(${PROJECT_NAME} main.cpp definitions.h)

set_common_target_properties(${PROJECT_NAME} "HERMES2D")

