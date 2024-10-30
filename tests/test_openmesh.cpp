#include "gtest/gtest.h"
#include "persistence.h"

TEST(OpenMeshTests, TestIfTriangleMesh) {
    const auto mesh = persistence::MeshDao::load_mesh_from_file("../tests/resources/fandisk.obj");

    for (auto it_f = mesh.faces_begin(); it_f != mesh.faces_end(); ++it_f) {
        EXPECT_EQ(3, mesh.valence(*it_f));
    }
}