cmake_minimum_required(VERSION 3.24)
project(MeshBound LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(Boost_USE_STATIC_LIBS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_LIBS)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")  # enable assert
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")  # enable assert
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS}")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG}")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-int-in-bool-context -Wno-sign-compare")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS}")
    set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")
endif ()

if (MSVC)
    # Enable M_PI and disable fopen() etc. warnings
    ADD_DEFINITIONS(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS)
    # Disable the following warnings : C4267 and 4244 (conversions) C4018 (sign compare) and C4800 (bool to int)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /wd4267 /wd4244 /wd4018 /wd4800")
endif ()

find_package(Eigen REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenMesh REQUIRED)

include(FetchContent)
include_directories(include)
include_directories(src)
include_directories(libs/pcg32/include)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN_INCLUDE_DIRS})

FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.14.1
)
FetchContent_MakeAvailable(spdlog)

# Targets
add_executable(
        meshbound

        src/dedge.cpp
        src/hierarchy.cpp
        src/main.cpp
        src/optimizer.cpp
        src/parametrizer.cpp
        src/persistence.cpp
        src/services.cpp
        src/bootstrap.cpp
        src/adapters.cpp
)

target_link_libraries(meshbound PRIVATE OpenMeshCore OpenMeshTools spdlog::spdlog)

# Testing
# For Windows: Prevent overriding the parent project's compiler/linker settings

FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
enable_testing()

add_executable(
        meshboundTests

        tests/test_persistence.cpp
        tests/test_openmesh.cpp
        src/dedge.cpp
        src/hierarchy.cpp
        src/optimizer.cpp
        src/parametrizer.cpp
        src/persistence.cpp
)


target_link_libraries(
        meshboundTests
        PRIVATE
        GTest::gtest_main
        OpenMeshCore
        OpenMeshTools
)

# Testing

include(GoogleTest)
gtest_discover_tests(meshboundTests)

