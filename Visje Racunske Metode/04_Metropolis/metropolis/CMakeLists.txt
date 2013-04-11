project(metropolis)

add_executable(metropolis main.cpp)
target_link_libraries(metropolis nutils)

install(TARGETS metropolis RUNTIME DESTINATION bin)
