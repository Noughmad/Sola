
include_directories(${CMAKE_CURRENT_BINARY_DIR})
qt4_automoc(integratortest.cpp)
add_executable(testgrate integratortest.cpp)
target_link_libraries(testgrate grate ${QT_QTCORE_LIBRARY} ${QT_QTTEST_LIBRARY})

add_test(testgrate testgrate)