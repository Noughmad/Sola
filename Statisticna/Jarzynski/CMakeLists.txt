project(jarzynski)

# OpenMP slows down execution, because rand() is no thread safe or reentrant
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(jarzynski system.cpp main.cpp)

install(TARGETS jarzynski RUNTIME DESTINATION bin)
