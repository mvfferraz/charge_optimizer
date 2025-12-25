cmake_minimum_required(VERSION 3.15)
project(ChargeOptimizer VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 3.3 QUIET)

if(NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found via find_package, will download automatically")
    
    include(FetchContent)
    FetchContent_Declare(
        Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
    )
    
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
    
    FetchContent_MakeAvailable(Eigen)
    
    add_library(Eigen3::Eigen ALIAS eigen)
else()
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
endif()

include_directories(${CMAKE_SOURCE_DIR}/src)

set(SOURCES
    src/main.cpp
    src/core/molecule.cpp
    src/core/esp_grid.cpp
    src/solver/qp_solver.cpp
    src/solver/active_set.cpp
    src/solver/constraints.cpp
    src/io/xyz_parser.cpp
    src/io/cube_parser.cpp
    src/analysis/validator.cpp
    src/analysis/symmetry.cpp
)

add_executable(charge_optimizer ${SOURCES})

target_link_libraries(charge_optimizer PRIVATE Eigen3::Eigen)

install(TARGETS charge_optimizer DESTINATION bin)

enable_testing()
add_subdirectory(tests)

message(STATUS "")
message(STATUS "Configuration Summary:")
message(STATUS "  Build type:        ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ compiler:      ${CMAKE_CXX_COMPILER}")
message(STATUS "  C++ flags:         ${CMAKE_CXX_FLAGS}")
message(STATUS "  Install prefix:    ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
