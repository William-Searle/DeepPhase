cmake_minimum_required(VERSION 3.11)
project(DeepPhase)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O2 -march=native -Wall -Wextra -Wpedantic)

# Find packages
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# Get NumPy include path
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NumPy_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Found NumPy at: ${NumPy_INCLUDE_DIR}")

# Gather sources
file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "include/*.h")

# Create a library target
add_library(DeepPhaseLib STATIC ${SOURCES} ${HEADERS})

target_include_directories(DeepPhaseLib PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${GSL_INCLUDE_DIRS}
    ${Python3_INCLUDE_DIRS}
    ${NumPy_INCLUDE_DIR}
)

target_link_libraries(DeepPhaseLib PUBLIC
    OpenMP::OpenMP_CXX
    GSL::gsl
    GSL::gslcblas
    Python3::Python
)

# Optional profiler linking
find_library(PROFILER_LIB profiler)
if(PROFILER_LIB)
    target_link_libraries(DeepPhaseLib PUBLIC ${PROFILER_LIB})
else()
    message(WARNING "Profiler library not found.")
endif()

# Create the standalone executable using the library
add_executable(deepphase_exe ${SOURCES})
target_link_libraries(deepphase_exe PRIVATE DeepPhaseLib)