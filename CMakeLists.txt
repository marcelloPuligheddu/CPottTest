cmake_minimum_required(VERSION 3.16)
project(CellularPotts)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(cpm
    src/CellularPotts.cpp
)

target_include_directories(cpm PUBLIC src)

# Fetch GoogleTest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(unit_test
    test/test_energy.cpp
    test/test_volume.cpp
)

target_link_libraries(unit_test
    cpm
    gtest_main
)

include(GoogleTest)
gtest_discover_tests(unit_test)
