cmake_minimum_required(VERSION 3.11)
project(DeepPhase VERSION 1.0.0 LANGUAGES CXX)

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O2 -march=native -Wall -Wextra -Wpedantic)

# Dependencies
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

# NumPy
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NumPy_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "Found NumPy at: ${NumPy_INCLUDE_DIR}")

# Source files
file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "include/*.h")

# Library target
add_library(DeepPhaseLib STATIC ${SOURCES} ${HEADERS})

target_include_directories(DeepPhaseLib PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
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

# Profiler (optional)
find_library(PROFILER_LIB profiler)
if(PROFILER_LIB)
    target_link_libraries(DeepPhaseLib PUBLIC ${PROFILER_LIB})
else()
    message(WARNING "Profiler library not found.")
endif()

# Install rules
include(GNUInstallDirs)

install(TARGETS DeepPhaseLib
    EXPORT DeepPhaseTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT DeepPhaseTargets
    FILE DeepPhaseTargets.cmake
    NAMESPACE DeepPhase::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DeepPhase
)

# Config + Version files
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/DeepPhaseConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/DeepPhaseConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/DeepPhaseConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DeepPhase
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/DeepPhaseConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/DeepPhaseConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/DeepPhase
)