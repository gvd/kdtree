project(kdtree)
cmake_minimum_required(VERSION 2.6)

list(APPEND CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")

find_package(Boost 1.50.0 REQUIRED COMPONENTS)
mark_as_advanced(Boost_DIR)

include_directories(SYSTEM ${Boost_INCLUDE_DIR})

add_executable(kdtree
    example.cpp
)
target_link_libraries(kdtree
)