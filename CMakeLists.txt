cmake_minimum_required(VERSION 3.1)
project(QuadriFlow LANGUAGES CXX)

option(BUILD_LOG "Enable verbose log" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(Boost_USE_STATIC_LIBS ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")  # enable assert
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")  # enable assert
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS}")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG}")

find_package(Eigen REQUIRED)
find_package(Boost REQUIRED)

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


include_directories(src)
include_directories(3rd/pcg32)
include_directories(3rd/pss)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${GLM_INCLUDE_DIRS})

if (BUILD_LOG)
    add_definitions(-DLOG_OUTPUT)
endif (BUILD_LOG)

set(
        quadriflow_SRC
        src/adjacent-matrix.cpp
        src/adjacent-matrix.hpp
        src/compare-key.hpp
        src/config.hpp
        src/dedge.cpp
        src/dedge.hpp
        src/disajoint-tree.hpp
        src/dset.hpp
        src/field-math.hpp
        src/flow.hpp
        src/hierarchy.cpp
        src/hierarchy.hpp
        src/loader.cpp
        src/loader.hpp
        src/localsat.cpp
        src/localsat.hpp
        src/main.cpp
        src/merge-vertex.cpp
        src/merge-vertex.hpp
        src/optimizer.cpp
        src/optimizer.hpp
        src/parametrizer.cpp
        src/parametrizer-flip.cpp
        src/parametrizer-int.cpp
        src/parametrizer-mesh.cpp
        src/parametrizer-scale.cpp
        src/parametrizer-sing.cpp
        src/parametrizer.hpp
        src/serialize.hpp
        src/subdivide.cpp
        src/subdivide.hpp
)

add_executable(
        quadriflow
        ${quadriflow_SRC}
)