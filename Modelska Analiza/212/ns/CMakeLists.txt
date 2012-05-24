project(ns)

add_executable(ns main.cpp)
target_link_libraries(ns gsl gslcblas)

