#include "gtest/gtest.h"
#include "spdlog/spdlog.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

#include "bootstrap.h"
#include "mathext.h"

using namespace Eigen;

std::vector<std::tuple<MatrixXf, MatrixXf, Vector3f, float> >
load_test_data() {
    std::ifstream inputFile("../tests/resources/svd-component.json");
    if (!inputFile.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
        throw std::runtime_error("Could not open the file!");
    }

    nlohmann::json jsonData;
    inputFile >> jsonData;

    std::vector<std::tuple<MatrixXf, MatrixXf, Vector3f, float> > examples;
    for (auto example: jsonData) {
        auto vertices = MatrixXf(example["vertices"].size(), 3);
        for (int i = 0; i < example["vertices"].size(); ++i) {
            const auto v = Vector3f(
                example["vertices"][i][0],
                example["vertices"][i][1],
                example["vertices"][i][2]
            );
            vertices.row(i) = v;
        }
        auto normals = MatrixXf(example["normals"].size(), 3);
        for (int i = 0; i < example["vertices"].size(); ++i) {
            const auto n = Vector3f(
                example["normals"][i][0],
                example["normals"][i][1],
                example["normals"][i][2]
            );
            normals.row(i) = n;
        }
        auto intersection = Vector3f(example["intersection"][0], example["intersection"][1],
                                     example["intersection"][2]);
        float error = example["error"];


        examples.emplace_back(std::make_tuple(
            vertices,
            normals,
            intersection,
            error
        ));
    }

    return examples;
}

TEST(SVDSuite, Case1) {
    auto vertices = MatrixXf(1, 3);
    vertices << 2, 0, 0;
    auto normals = MatrixXf(1, 3);
    normals << 1, 0, 0;

    const Vector4f result = mathext::intersect_planes(vertices, normals);

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

    const Vector4f result = mathext::intersect_planes(vertices, normals);

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

    const Vector4f result = mathext::intersect_planes(vertices, normals);

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

    const Vector4f result = mathext::intersect_planes(vertices, normals);

    const auto expected = Vector4f(3, 0, 0, 0);
    ASSERT_TRUE(result.isApprox(expected, 1e-2));
}

TEST(SVDSuite, ComparisonExamples) {
    auto examples = load_test_data();

    int i = 0;
    auto inputs = MatrixXf(examples.size(), 6);
    auto result = MatrixXf(examples.size(), 8);
    for (auto &[vertices, normals, x, error]: examples) {
        Vector3f v_sum = vertices.colwise().sum();
        Vector3f n_sum = normals.colwise().sum();
        inputs.block<1, 3>(0, 0) = v_sum;
        inputs.block<1, 3>(0, 3) = n_sum;

        const auto xerr = mathext::intersect_planes(vertices, normals);
        result.row(i).block<1, 3>(0, 0) = x;
        result.row(i).block<1, 3>(0, 3) = xerr.head<3>();
        result.row(i)[6] = error;
        result.row(i)[7] = xerr[3];
        i++;
    }

    std::cout << "Inputs:\n" << inputs.block<10, 6>(0, 0) << std::endl;
    std::cout << "Results:\n" << result.block<10, 6>(0, 0) << std::endl;
    std::cout << "Difference:\n" << result.block<10, 3>(0, 0) - result.block<10, 3>(0, 3) << std::endl;

    ASSERT_TRUE(result.block(0, 0, result.rows(), 3).isApprox(result.block(0, 3, result.rows(), 3), 1e-1));
}

TEST(MathExtensions, NormalizeAndDenormalize) {
    const auto mesh = bootstrap::Container().mesh_service().load_mesh("../tests/resources/box.ply");

    MatrixXd vertices(3, mesh.n_vertices());
    for (int i = 0; i < mesh.n_vertices(); ++i) {
        auto v = mesh.point(entities::Mesh::VertexHandle(i));
        vertices.col(i) = Vector3d(v[0], v[1], v[2]);
    }

    const auto [vertices_norm, scale, offset] = mathext::normalize(vertices);
    const MatrixXd vertices_back = mathext::denormalize(vertices_norm, scale, offset);

    EXPECT_TRUE(vertices.isApprox(vertices_back, 1e-6));
}
