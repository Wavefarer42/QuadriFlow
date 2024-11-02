#include "gtest/gtest.h"
#include "persistence.h"

const std::string path_mesh = "../tests/resources/fandisk.obj";
const std::string path_model = "../tests/resources/Bear_2.ubs";

TEST(MeshDaoTests, LoadMeshFromFileObj) {
    const persistence::MeshDao sut;

    const auto result = sut.load_mesh_from_file(path_mesh);

    EXPECT_EQ(7229, result.n_vertices());
    EXPECT_EQ(14454, result.n_faces());
}


TEST(MeshDaoTests, LoadUnboundSdfn) {
    const persistence::MeshDao sut;

    const auto result = sut.load_sdfn_from_file(path_model);

    EXPECT_EQ(50, result.size());
}

TEST(MeshDaoTests, SampleUnboundSdfn) {
    const persistence::MeshDao dao;

    const auto sut = dao.load_sdfn_from_file(path_model)[0];

    Eigen::MatrixXf domain(5, 3);
    domain << 1, 0, 0,
            0, 1, 0,
            0, 1, 1,
            1, 1, 0,
            0, 0, 1;

    const Eigen::MatrixXf result = sut(domain);
    EXPECT_NEAR(result.norm(), 11.6942, 1e-4);
}