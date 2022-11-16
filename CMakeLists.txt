# =============================================================================
# BDcpp -- Simple Bjontegaard Delta metric implementation for C++.
# 
# MIT License
#
# Copyright (c) 2022 Tim Bruylants
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# =============================================================================

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(bdcpp LANGUAGES CXX C)

# Setup Conan package finding.
set(CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})

# Language settings.
set(CMAKE_CXX_STANDARD 17)  # C++17
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# # Some global settings.
# if(MSVC)
#   # Disable some warnings on MSVC.
#   add_compile_definitions(_SCL_SECURE_NO_WARNINGS)
#   add_compile_options("/wd\"26451\"")
# endif()

if (MSVC)
    add_compile_options(/W4 /bigobj "/wd\"4127\"")
endif()

# Get everything into one final folder.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# External packages.
find_package(Eigen3 CONFIG)

# -----------------------------------
# Library code.
# -----------------------------------
set(SRC_BDCPP_LIB_PRIVATE
    bdcpp/bdcpp.cpp
    bdcpp/bdcpp.h
)

add_library(bdcpp STATIC ${SRC_BDCPP_LIB_PRIVATE})
target_link_libraries(bdcpp PRIVATE Eigen3::Eigen)
target_include_directories(bdcpp PUBLIC "bdcpp/")

# -----------------------------------
# Main example application.
# -----------------------------------
set(SRC_FILES
    example.cpp
)
add_executable(example ${SRC_FILES})
target_link_libraries(example PRIVATE bdcpp)
