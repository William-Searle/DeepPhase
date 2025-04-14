cmake_minimum_required(VERSION 3.11)
project(ISSMGW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create executable
include_directories(src)
file(GLOB SOURCES "src/*.cpp")
add_executable(deepphase ${SOURCES})

message(STATUS "Sources: ${SOURCES}")

# setting warning compiler flags
add_compile_options(-Wall -Wextra -Wpedantic)

# specify install rules
install(TARGETS deepphase DESTINATION bin)