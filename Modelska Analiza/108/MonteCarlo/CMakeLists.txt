project(montecarlo)

set(CMAKE_BUILD_TYPE Release)

add_executable(montecarlo main.cpp)
target_link_libraries(montecarlo gsl gslcblas)

