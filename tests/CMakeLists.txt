include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
link_libraries(conf)

# Add all test source files which should be compiled into executables,
# functions below add the rest.
set (test_source_files
     test_system_conf.cpp
     test_md_integrator.cpp)

# Files which are used for tests need the correct path inside of the
# tests. The base path is set in this variable during configuration.
set (Test_Files_Directory ${PROJECT_BINARY_DIR}/tests/include)
configure_file (
    "${PROJECT_SOURCE_DIR}/tests/utils.h.in"
    "${PROJECT_BINARY_DIR}/tests/utils.h"
)

# Find and copy the required test files to the correct path.
file(GLOB test_files ${PROJECT_SOURCE_DIR}/tests/include/*.gro)
message(${test_files})
add_custom_target(copy_test_files)
foreach(file ${test_files})
    add_custom_command(TARGET copy_test_files POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E
                           make_directory ${Test_Files_Directory})
    add_custom_command(TARGET copy_test_files POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E
                           copy ${file} ${Test_Files_Directory})
endforeach()

foreach (file IN ITEMS ${test_source_files})
    get_filename_component(target ${file} NAME_WE)
    add_executable(${target} ${file})
    add_dependencies(${target} copy_test_files)
    add_test(NAME ${target} COMMAND ${target})
endforeach(file)
