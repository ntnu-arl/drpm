cmake_minimum_required(VERSION 3.0.2)
project(drpm)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED NO_MODULE) 

# Get Sophus (TODO: Remove if not used)
include(FetchContent)
set(SOPHUS_USE_BASIC_LOGGING ON CACHE BOOL "Don't use fmt for Sophus")
set(BUILD_SOPHUS_TESTS OFF CACHE BOOL "Don't build Sophus tests")
set(BUILD_SOPHUS_EXAMPLES OFF CACHE BOOL "Don't build Sophus Examples")
FetchContent_Declare(sophus SYSTEM URL https://github.com/strasdat/Sophus/archive/refs/tags/1.22.10.tar.gz)
FetchContent_MakeAvailable(sophus)


add_executable(example
            src/example.cpp)
target_include_directories(example PRIVATE ${PROJECT_SOURCE_DIR}/src ${Boost_INCLUDE_DIRS})

target_link_libraries(example Eigen3::Eigen Sophus::Sophus)