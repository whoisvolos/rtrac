cmake_minimum_required(VERSION 3.0.2)
project(rtrac)

set(CMAKE_INCLUDE_PATH "")

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# Find the QtWidgets library
set(CMAKE_PREFIX_PATH /usr/local/opt/qt5)

find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)
find_package(OpenGL REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BUILD_TYPE Release)

set(SOURCE_FILES
    main.cpp
    funcs.hpp
    time.hpp
    integrate.hpp
    integrate.cpp
    math/operations.h
    math/operations_ext.h
    math/plane.h
    math/plane.cpp
    math/triangle.h
    math/triangle.cpp
    math/types.h
    gui/main_window.cpp
    gui/main_window.h
    )

qt5_wrap_ui(UI_HEADERS gui/main_window.ui)

message(STATUS ${UI_HEADERS})
add_executable(rtrac ${SOURCE_FILES} ${UI_HEADERS})

qt5_use_modules(rtrac Core Gui Widgets OpenGL)

target_link_libraries(rtrac Qt5::Widgets Qt5::OpenGL ${OPENGL_LIBRARIES})

# copy ui_* for IDE intellisense
ADD_CUSTOM_COMMAND(TARGET rtrac
          POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_BINARY_DIR}/ui_main_window.h" ${PROJECT_SOURCE_DIR}
)
