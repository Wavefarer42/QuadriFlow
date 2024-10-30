cmake_minimum_required(VERSION 3.24)

add_executable(
        meshbound

        test_migration.cpp
)


target_link_libraries(meshbound GTest::gtest_main)

# Testing

include(GoogleTest)
gtest_discover_tests(meshbound)