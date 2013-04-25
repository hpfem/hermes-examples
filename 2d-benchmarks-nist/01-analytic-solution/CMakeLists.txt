project(nist-01) 
add_executable(${PROJECT_NAME} main.cpp ../NIST-util.cpp definitions.cpp)
set_common_target_properties(${PROJECT_NAME} "HERMES2D")  


