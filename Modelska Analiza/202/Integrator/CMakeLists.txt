project(integrator)
find_package(Qt4 COMPONENTS QtCore QtGui REQUIRED)
cmake_minimum_required(VERSION 2.8)

enable_testing()

include_directories({${QT_INCLUDES})

set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++0x)

add_executable(integrator main.cpp)
target_link_libraries(integrator grate ${QT_QTCORE_LIBRARY})

add_subdirectory(lib)

add_subdirectory(tests)