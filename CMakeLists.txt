project(kdtree)
cmake_minimum_required(VERSION 2.6)

list(APPEND CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")

find_package(Boost 1.50.0 REQUIRED COMPONENTS)
mark_as_advanced(Boost_DIR)

include_directories(SYSTEM ${Boost_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

find_package(GTest)
if(GTEST_FOUND)
    include(FindThreads)
    if(Threads_FOUND)
        set(GTEST_LIBRARIES ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
        enable_testing()
        add_subdirectory(test)
    endif(Threads_FOUND)
else(GTEST_FOUND)
    message(WARNING "GTest not found, not compiling tests")
endif(GTEST_FOUND)