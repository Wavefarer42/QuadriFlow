#include "gtest/gtest.h"
#include "persistence.h"

TEST(OpenMeshTests, TestIfTriangleMesh) {
    persistence::MeshDao sut;
    const auto mesh = sut.load_mesh_from_file("../tests/resources/fandisk.obj");

    for (auto it_f = mesh.faces_begin(); it_f != mesh.faces_end(); ++it_f) {
        EXPECT_EQ(3, mesh.valence(*it_f));
    }
}