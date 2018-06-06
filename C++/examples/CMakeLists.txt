# PROJECT CONFIGURATION
cmake_minimum_required(VERSION 3.1)

# Set this project's configuration
project(SESyncExecutable C CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON) # We require C++ 14

find_library(SESync sesync REQUIRED)

# SE-Sync command-line driver
add_executable(SE-Sync main.cpp)
target_link_libraries(SE-Sync SESync)
if(${ENABLE_PROFILING})
target_compile_definitions(SE-Sync PRIVATE GPERFTOOLS)
target_link_libraries(SE-Sync profiler)
endif()


message(STATUS "Building main SE-Sync command-line executable in directory ${EXECUTABLE_OUTPUT_PATH}\n")
