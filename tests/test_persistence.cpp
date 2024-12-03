#include "gtest/gtest.h"
#include "persistence.h"
#include "bootstrap.h"

const std::string path_mesh = "../tests/resources/fandisk.obj";
const std::string path_model = "../tests/resources/Bear_2.ubs";

TEST(MeshDaoTests, LoadMeshFromFileObj) {
    const persistence::MeshDao sut = bootstrap::Container().mesh_dao();

    const auto result = sut.load_mesh(path_mesh);

    EXPECT_EQ(7229, result.n_vertices());
    EXPECT_EQ(14454, result.n_faces());
}


TEST(MeshDaoTests, LoadUnboundSdfn) {
    const persistence::MeshDao sut = bootstrap::Container().mesh_dao();

    const auto result = sut.load_model(path_model);

    EXPECT_EQ(50, result.size());
}

