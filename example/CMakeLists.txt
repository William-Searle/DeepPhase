# --------------------------
# CMakeLists.txt for DeepPhase/example
# --------------------------

# -------------------------
# Output directory for binaries
# -------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# --------------------------
# Link libraries
# --------------------------
link_libraries(
    DeepPhaseLib
)

# --------------------------
# Example executables
# --------------------------
add_executable(run_profile run_profile.cpp)