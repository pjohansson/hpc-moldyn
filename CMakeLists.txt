cmake_minimum_required (VERSION 3.0)
project ("HPC Moldyn")

enable_testing()
add_subdirectory(tests)

add_compile_options(-Wall -std=c++14 -O3)

add_executable(md src/main.cpp)

add_library(analytics STATIC src/analytics.cpp)
add_library(conf STATIC src/conf.cpp)
add_library(integrator STATIC src/integrator.cpp)
add_library(mpi STATIC src/mpi_impl.cpp)
add_library(params STATIC src/params.cpp)

target_link_libraries(md analytics conf integrator mpi params)
