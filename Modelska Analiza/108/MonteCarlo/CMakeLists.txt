project(montecarlo)

add_executable(montecarlo main.cpp)
target_link_libraries(montecarlo gsl gslcblas)

