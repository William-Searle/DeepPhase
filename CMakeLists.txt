cmake_minimum_required(VERSION 3.11)
project(ISSMGW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -march=native")

# Add ALGLIB static library
file(GLOB ALGLIB_SRC "include/alglib/*.cpp")
if(NOT ALGLIB_SRC)
    message(FATAL_ERROR "No source files found in alglib/. Check your path or GLOB pattern.")
endif()
add_library(alglib STATIC ${ALGLIB_SRC})
target_include_directories(alglib PUBLIC ${CMAKE_SOURCE_DIR}/include/alglib)

# omp parallelisation
find_package(OpenMP REQUIRED)

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
file(GLOB SOURCES "src/*.cpp")
add_executable(deepphase ${SOURCES})

target_include_directories(deepphase PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)

message(STATUS "Sources: ${SOURCES}")

# profiler - add dev version so it doesn't always run
find_library(PROFILER_LIB profiler)

target_link_libraries(deepphase PRIVATE 
    OpenMP::OpenMP_CXX
    alglib
    ${Python3_LIBRARIES}
    ${PROFILER_LIB}
)
message(STATUS "Successfully linked to dependent libraries")

# setting warning compiler flags
add_compile_options(-Wall -Wextra -Wpedantic)

# specify install rules
install(TARGETS deepphase DESTINATION bin)