project(dft)

add_executable(dft main.cpp)
target_link_libraries(dft gsl gslcblas)

install(TARGETS dft RUNTIME DESTINATION bin)
