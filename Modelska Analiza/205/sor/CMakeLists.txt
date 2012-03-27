project(sor)

add_executable(sor main.cpp)
target_link_libraries(sor gsl gslcblas)

