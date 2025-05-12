cmake_minimum_required(VERSION 3.11)
project(ISSMGW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

# Find GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

# Find Python
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
include_directories(${Python3_INCLUDE_DIRS})

# Find numpy
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NumPy_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
include_directories(${NumPy_INCLUDE_DIR})
message(STATUS "Found NumPy at: ${NumPy_INCLUDE_DIR}")

# Create executable
include_directories(src)
file(GLOB SOURCES "src/*.cpp")
add_executable(deepphase ${SOURCES})

message(STATUS "Sources: ${SOURCES}")

# profiler - add dev version so it doesn't always run
find_library(PROFILER_LIB profiler)

target_link_libraries(deepphase PRIVATE 
    gsl 
    gslcblas 
    m
    ${Python3_LIBRARIES}
    ${PROFILER_LIB}
)
message(STATUS "Successfully linked to dependent libraries")

# setting warning compiler flags
add_compile_options(-Wall -Wextra -Wpedantic)

# specify install rules
install(TARGETS deepphase DESTINATION bin)