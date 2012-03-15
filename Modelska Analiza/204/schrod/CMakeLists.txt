project(schrod)

add_executable(schrod main.cpp)
target_link_libraries(schrod gsl gslcblas)

