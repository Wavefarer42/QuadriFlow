cmake_minimum_required(VERSION 3.24)
project(MeshBound LANGUAGES CXX)

option(ENABLE_DEV_DEBUG "Enable development debugging artifacts" ON)

if(ENABLE_DEV_DEBUG)
    add_compile_definitions(DEV_DEBUG=1)
else()
    add_compile_definitions(DEV_DEBUG=0)
endif()

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


include(FetchContent)
include_directories(include)
include_directories(src)

# Tbb
find_package(TBB REQUIRED)

# Eigen
find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})

# Boost
find_package(Boost REQUIRED COMPONENTS thread)
include_directories(${Boost_INCLUDE_DIRS})

# OpenMesh
find_package(OpenMesh REQUIRED)

# PCG32
include_directories(third_party/pcg32/include)

# Spdlog
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.14.1
)
FetchContent_MakeAvailable(spdlog)
include_directories(${spdlog_SOURCE_DIR})

# GLM
FetchContent_Declare(
        glm
        GIT_REPOSITORY https://github.com/g-truc/glm.git
        GIT_TAG 0.9.8  # Use a specific version here (latest tag is preferable)
)
FetchContent_MakeAvailable(glm)
include_directories(${glm_SOURCE_DIR})

# Nlohmann JSON
FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.2  # Specify the version you want to use
)
FetchContent_MakeAvailable(nlohmann_json)
include_directories(${nlohmann_json_SOURCE_DIR})

# Argparse
include(FetchContent)
FetchContent_Declare(
        argparse
        GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

# ----- UnboundEval -----

# Unbound
set(UNBOUNDEVAL_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/UnboundEval/include)
set(UNBOUNDEVAL_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/UnboundEval/src)

include_directories(third_party/UnboundEval/include)
include_directories(third_party/UnboundEval/include/UB)
include_directories(third_party/UnboundEval/include/Core)
include_directories(third_party/UnboundEval/include/xxhash)
include_directories(third_party/UnboundEval/include/stduuid)
include_directories(third_party/UnboundEval/include/Math)

add_library(
        UnboundEval STATIC
        ${UNBOUNDEVAL_SRC_DIR}/Collection.cpp
        ${UNBOUNDEVAL_SRC_DIR}/Collection.cpp
        ${UNBOUNDEVAL_SRC_DIR}/Edit.cpp
        ${UNBOUNDEVAL_SRC_DIR}/Instruction.cpp
        ${UNBOUNDEVAL_SRC_DIR}/Shape.cpp
        ${UNBOUNDEVAL_SRC_DIR}/UB.cpp
        ${UNBOUNDEVAL_SRC_DIR}/Centripetal.cpp
        third_party/UnboundEval/src/xxhash.cpp
)

target_compile_definitions(
        UnboundEval
        PRIVATE
        GLM_ENABLE_EXPERIMENTAL
)

# ----- MeshBound ----

# Targets
add_executable(
        meshbound

        src/main.cpp

        src/dedge.cpp
        src/hierarchy.cpp
        src/optimizer.cpp
        src/parametrizer.cpp
        src/persistence.cpp
        src/services.cpp
        src/bootstrap.cpp
        src/adapters.cpp
        src/surfacenets.cpp
)

target_link_libraries(
        meshbound
        PRIVATE
        OpenMeshCore
        OpenMeshTools
        spdlog::spdlog
        nlohmann_json::nlohmann_json
        UnboundEval
        "-framework CoreFoundation"
        TBB::tbb
        argparse
)

# Testing
# For Windows: Prevent overriding the parent project's compiler/linker settings

FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
enable_testing()

# Build Targets
add_executable(
        meshboundTests

        tests/test_persistence.cpp
        tests/test_libraries.cpp
        tests/test_surfacenets.cpp

        src/dedge.cpp
        src/hierarchy.cpp
        src/optimizer.cpp
        src/parametrizer.cpp
        src/persistence.cpp
        src/services.cpp
        src/bootstrap.cpp
        src/adapters.cpp
        src/surfacenets.cpp
)

target_link_libraries(
        meshboundTests
        PRIVATE
        GTest::gtest_main
        OpenMeshCore
        OpenMeshTools
        spdlog::spdlog
        nlohmann_json::nlohmann_json
        UnboundEval
        "-framework CoreFoundation"
        TBB::tbb
)

# Testing

include(GoogleTest)
gtest_discover_tests(meshboundTests)

