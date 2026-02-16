cmake_minimum_required(VERSION 3.12.0)
project(nmpde-wave-equation LANGUAGES CXX C)

include(cmake-common.cmake)

add_library(WaveEquation STATIC src/WaveEquation.cpp)
deal_ii_setup_target(WaveEquation)

add_executable(wave-equation src/main.cpp)
target_link_libraries(wave-equation WaveEquation)
deal_ii_setup_target(wave-equation)

if (APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_include_directories(WaveEquation BEFORE PRIVATE
        "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1")
    target_include_directories(wave-equation BEFORE PRIVATE
        "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1")
endif ()
