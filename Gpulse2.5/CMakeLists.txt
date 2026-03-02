cmake_minimum_required(VERSION 3.10)
project(Gpulse VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O3 -march=native -ffast-math)

# Find OpenMP (Required for parallel processing)
find_package(OpenMP REQUIRED)

# Find HDF5 (Required for export tool only)
find_package(HDF5 COMPONENTS CXX)

# Include directories (Public headers)
include_directories(include)

# Define the shared library
add_library(gpulse_lib SHARED
    src/Gpulse.cpp
    src/constants.cpp
)

# Link OpenMP to the library
target_link_libraries(gpulse_lib PUBLIC OpenMP::OpenMP_CXX)

# Define the test executable
add_executable(test_gpulse tests/test_gpulse.cpp)
target_link_libraries(test_gpulse PRIVATE gpulse_lib)
set_target_properties(test_gpulse PROPERTIES BUILD_RPATH "\${CMAKE_BINARY_DIR}")

# Define the export tool (only if HDF5 is found)
if(HDF5_FOUND)
    add_executable(export_field tools/export_field.cpp)
    target_link_libraries(export_field PRIVATE gpulse_lib HDF5::HDF5)
    set_target_properties(export_field PROPERTIES BUILD_RPATH "\${CMAKE_BINARY_DIR}")
    
    # Copy plot_field.py to build directory
    configure_file(tools/plot_field.py ${CMAKE_BINARY_DIR}/plot_field.py COPYONLY)

    # Copy plot_temporal.py to build directory
    configure_file(tools/plot_temporal.py ${CMAKE_BINARY_DIR}/plot_temporal.py COPYONLY)
    
    # Temporal Power Exporter Tool
    add_executable(export_temporal tools/export_temporal.cpp)
    target_link_libraries(export_temporal PRIVATE gpulse_lib HDF5::HDF5)
    set_target_properties(export_temporal PROPERTIES BUILD_RPATH "\${CMAKE_BINARY_DIR}")

    message(STATUS "HDF5 found. 'export_field' tool will be built.")
else()
    message(WARNING "HDF5 not found. 'export_field' tool will NOT be built.")
endif()
