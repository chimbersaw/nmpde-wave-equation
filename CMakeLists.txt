cmake_minimum_required(VERSION 3.12.0)
project(nmpde-wave-equation LANGUAGES CXX C)

include(cmake-common.cmake)

add_library(Poisson1D STATIC src/Poisson1D.cpp)
deal_ii_setup_target(Poisson1D)

add_executable(wave-equation src/main.cpp)
target_link_libraries(wave-equation Poisson1D)
deal_ii_setup_target(wave-equation)

if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_include_directories(Poisson1D BEFORE PRIVATE
        "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1")
    target_include_directories(wave-equation BEFORE PRIVATE
        "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1")
endif()
