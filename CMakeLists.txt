cmake_minimum_required(VERSION 3.20)
project(esercizi_lsn CXX)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
#set(CMAKE_OUT_DIR "")


add_library(project_warnings INTERFACE)
add_library(project_options INTERFACE)
target_compile_options(project_warnings
        INTERFACE
        -Wall
        -Wextra
        -Wshadow
        -Wnon-virtual-dtor
        -Wold-style-cast
        -Wcast-align
        -Woverloaded-virtual
        -Wpedantic
        -Wconversion
        -Wno-sign-conversion
        -Wnull-dereference
        -Wdouble-promotion
        -Wformat=2
        -Werror
        -Wno-unused-variable
        -Wno-psabi)


add_subdirectory(src)
add_subdirectory(tests)
