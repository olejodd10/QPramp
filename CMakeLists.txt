cmake_minimum_required(VERSION 3.10)

project(QPramp VERSION 1.0)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${QPramp_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${QPramp_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QPramp_SOURCE_DIR})

add_library(compiler_flags INTERFACE)
target_compile_options(compiler_flags INTERFACE -Ofast -g)

add_subdirectory(src)
