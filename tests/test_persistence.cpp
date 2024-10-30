#include "gtest/gtest.h"
#include "persistence.h"

TEST(MeshDaoTests, LoadMeshFromFileObj) {
    const auto result = persistence::MeshDao::load_mesh_from_file("../tests/resources/fandisk.obj");
    EXPECT_EQ(7229, result.n_vertices());
    EXPECT_EQ(14454, result.n_faces());
}
