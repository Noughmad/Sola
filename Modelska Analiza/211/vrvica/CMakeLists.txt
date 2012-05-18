project(vrvica)

add_executable(vrvica main.cpp)
target_link_libraries(vrvica gsl gslcblas QtCore QtGui)

