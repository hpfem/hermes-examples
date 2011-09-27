project(wall-on-fire-adapt-space-and-time)
add_executable(${PROJECT_NAME} main.cpp definitions.cpp definitions.h)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")