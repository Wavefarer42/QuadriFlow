#include "gtest/gtest.h"
#include "bootstrap.h"


TEST(LibraryTests, QuadMeshLoading) {
    const persistence::MeshDao sut = bootstrap::Container().mesh_dao();

    const auto mesh = sut.load_mesh_from_file("../tests/resources/fandisk.obj");

    for (auto it_f = mesh.faces_begin(); it_f != mesh.faces_end(); ++it_f) {
        EXPECT_EQ(3, mesh.valence(*it_f));
    }
}

TEST(LibraryTests, EigenBroadcasting) {
    Eigen::MatrixXf lhs(4, 3); // Example: 4 rows, 3 columns
    lhs << 1, 2, 3,
            4, 5, 6,
            7, 8, 9,
            10, 11, 12;

    Eigen::VectorXf rhs(4); // Example vector with 4 elements (one per row)
    rhs << 1, 2, 3, 4;

    // Divide each row of the matrix by the corresponding entry in the vector
    Eigen::MatrixXf result = lhs.array().colwise() / rhs.array();
    Eigen::MatrixXf expected = Eigen::MatrixXf(4, 3);
    expected << 1, 2, 3,
            2, 2.5, 3,
            2.33, 2.66, 3,
            2.5, 2.75, 3;

    std::cout << "Result:\n" << result << std::endl;

    EXPECT_TRUE(result.isApprox(expected, 1e-2));
}

TEST(LibraryTests, EigenBroadcastingMultiplication) {
    Eigen::MatrixXf lhs(4, 3); // Example: 4 rows, 3 columns
    lhs << 1, 2, 3,
            1, 2, 3,
            1, 2, 3,
            1, 2, 3;

    Eigen::VectorXf rhs(4); // Example vector with 4 elements (one per row)
    rhs << 1, 2, 3, 4;

    // Divide each row of the matrix by the corresponding entry in the vector
    Eigen::MatrixXf result = lhs.array().colwise() * rhs.array();
    Eigen::MatrixXf expected = Eigen::MatrixXf(4, 3);
    expected << 1, 2, 3,
            2, 4, 6,
            3, 6, 9,
            4, 8, 12;

    std::cout << "Result:\n" << result << std::endl;

    EXPECT_TRUE(result.isApprox(expected, 1e-2));
}
