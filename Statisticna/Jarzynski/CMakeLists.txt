project(jarzynski)

add_executable(jarzynski system.cpp main.cpp)

install(TARGETS jarzynski RUNTIME DESTINATION bin)
