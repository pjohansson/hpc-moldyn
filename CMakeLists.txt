cmake_minimum_required (VERSION 3.0)
project ("HPC moldyn")

add_compile_options(-Wall -std=c++14 -O3)

add_executable(md src/main.cpp)

add_library(conf STATIC src/conf.cpp)
target_link_libraries(md conf)

enable_testing()
add_subdirectory(tests)
