cmake_minimum_required(VERSION 3.16)
project(fluidos_lab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Same as before: the binary ends up in the project root (./f1)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

option(USE_OPENMP "Parallelise the loops with OpenMP" ON)

# ---- raylib via pkg-config ----
find_package(PkgConfig REQUIRED)
pkg_check_modules(RAYLIB REQUIRED IMPORTED_TARGET raylib)

# ---- sources ----
add_executable(f1
    src/main.cpp
    src/BC.cpp
    src/cylinder.cpp
    src/mouse.cpp
    src/solver.cpp
    src/vorticity.cpp
    src/multigrid.cpp
)

target_include_directories(f1 PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(f1 PRIVATE PkgConfig::RAYLIB)
target_compile_options(f1 PRIVATE -Wall -Wextra -Wshadow)

# ---- OpenMP (optional) ----
if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(f1 PRIVATE OpenMP::OpenMP_CXX)
        message(STATUS "OpenMP: enabled")
    else()
        message(WARNING "OpenMP requested but not found; building single-threaded")
    endif()
else()
    # the #pragma omp lines are harmless but would warn with -Wall
    target_compile_options(f1 PRIVATE -Wno-unknown-pragmas)
endif()

# ---- tune for this machine (native build only) ----
if(CMAKE_BUILD_TYPE STREQUAL "Release" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(f1 PRIVATE -march=native)
endif()
