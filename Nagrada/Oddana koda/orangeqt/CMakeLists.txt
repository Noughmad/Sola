project(orangeqt)

cmake_minimum_required(VERSION 2.6)
find_package(Qt4 COMPONENTS QtCore QtGui REQUIRED)
find_package(PythonInterp REQUIRED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/modules)
find_package(PyQt4 REQUIRED)
find_package(SIP REQUIRED)
find_package(PythonLibrary REQUIRED)
include(SIPMacros)

include_directories(${QT_INCLUDES} ${SIP_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${PYTHON_INCLUDE_PATH})

#### orangeqt sources ####
set(orangeqt_SRCS 
    plotitem.cpp 
    point.cpp 
    curve.cpp 
    unconnectedlinescurve.cpp 
    networkcurve.cpp
    multicurve.cpp
    plot.cpp 
)
qt4_automoc(${orangeqt_SRCS})

#### Python bindings ####
set(SIP_PLATFORM_TAG)
if (Q_WS_WIN)
	set(PYQT4_PLATFORM_TAG "WS_WIN")
endif (Q_WS_WIN)
if (Q_WS_X11)
	set(PYQT4_PLATFORM_TAG "WS_X11")
endif (Q_WS_X11)
if (Q_WS_MAC)
	set(PYQT4_PLATFORM_TAG "WS_MACX")
endif (Q_WS_MAC)
set(SIP_TAGS ${PYQT4_VERSION_TAG} ${PYQT4_PLATFORM_TAG})
set(SIP_INCLUDES ${PYQT4_SIP_DIR})
set(SIP_CONCAT_PARTS 1)
set(SIP_EXTRA_OPTIONS "-k")

set(SIP_EXTRA_LINK_LIBRARIES ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTOPENGL_LIBRARY})
set(SIP_EXTRA_SOURCES ${orangeqt_SRCS})

add_sip_python_module(orangeqt orangeqt.sip)

#### copying the library after building ####
if (NOT ORANGE_LIB_DIR)
    set (ORANGE_LIB_DIR "../..")
endif()

get_target_property(PYTHON_LIB_LOCATION python_module_orangeqt LOCATION)

add_custom_command(
   TARGET python_module_orangeqt
   POST_BUILD
   COMMAND ${CMAKE_COMMAND} -E copy ${PYTHON_LIB_LOCATION} ${ORANGE_LIB_DIR}
 )
