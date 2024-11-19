#include "gtest/gtest.h"
#include "spdlog/spdlog.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <Eigen/Core>

#include "transformations.h"

using namespace Eigen;

TEST(SVDSuite, Case1) {
    auto vertices = MatrixXf(1, 3);
    vertices << 2, 0, 0;
    auto normals = MatrixXf(1, 3);
    normals << 1, 0, 0;

    const Vector4f result = transformations::intersect_planes(vertices, normals);

    const auto expected = Vector4f(2, 0, 0, 0);
    ASSERT_TRUE(result.isApprox(expected, 1e-2));
}

TEST(SVDSuite, Case2) {
    auto vertices = MatrixXf(2, 3);
    vertices << 1, 0, 0,
            0, 2, 0;
    auto normals = MatrixXf(2, 3);
    normals << 1, 0, 0,
            0, 1, 0;

    const Vector4f result = transformations::intersect_planes(vertices, normals);

    const auto expected = Vector4f(1, 2, 0, 0);
    ASSERT_TRUE(result.isApprox(expected, 1e-2));
}

TEST(SVDSuite, Case3) {
    auto vertices = MatrixXf(3, 3);
    vertices << 3, 0, 0,
            0, 4, 0,
            0, 0, 5;
    auto normals = MatrixXf(3, 3);
    normals << 1, 1, 0,
            0, 1, 1,
            1, 0, 1;

    const Vector4f result = transformations::intersect_planes(vertices, normals);

    const auto expected = Vector4f(2, 1, 3, 0);
    ASSERT_TRUE(result.isApprox(expected, 1e-2));
}

TEST(SVDSuite, Case4) {
    auto vertices = MatrixXf(3, 3);
    vertices << 1, 0, 0,
            3, 0, 0,
            5, 0, 0;
    auto normals = MatrixXf(3, 3);
    normals << 1, 0, 0,
            1, 0, 0,
            1, 0, 0;

    const Vector4f result = transformations::intersect_planes(vertices, normals);

    const auto expected = Vector4f(3, 0, 0, 0);
    ASSERT_TRUE(result.isApprox(expected, 1e-2));
}