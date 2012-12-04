project(kdtree_tests)

add_executable(kdtree_tests
    tests.cpp
)
target_link_libraries(kdtree_tests
    ${GTEST_LIBRARIES}
)
add_test(kdtree_tests kdtree_tests)
