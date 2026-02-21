cmake_minimum_required(VERSION 3.12.0)
project(nmpde-wave-equation LANGUAGES CXX C)

include(cmake-common.cmake)

add_executable(wave-equation
    src/main.cpp
    src/config.cpp
    src/function_factory.cpp
    src/wave_solver.cpp
    src/integrators/theta_integrator.cpp
    src/integrators/newmark_integrator.cpp
    src/convergence_runner.cpp)

deal_ii_setup_target(wave-equation)

if (APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_include_directories(wave-equation BEFORE PRIVATE
        "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1")
endif ()
