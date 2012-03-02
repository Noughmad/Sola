project(integrator)
find_package(Qt4 COMPONENTS QtCore REQUIRED)
cmake_minimum_required(VERSION 2.8)

enable_testing()

include_directories({${QT_INCLUDES})

add_executable(integrator main.cpp)
target_link_libraries(integrator grate ${QT_QTCORE_LIBRARY})

add_subdirectory(lib)

add_subdirectory(tests)