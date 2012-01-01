project(populacije)

set (CMAKE_CXX_FLAGS "-Wall -O3 -std=c++0x")

add_executable(populacije main.cpp)
target_link_libraries(populacije gsl gslcblas)
