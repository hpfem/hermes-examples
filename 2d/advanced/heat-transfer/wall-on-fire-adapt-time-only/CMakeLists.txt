project(wall-on-fire-adapt-time-only)
add_executable(${PROJECT_NAME} main.cpp definitions.cpp definitions.h)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")