#include "gtest/gtest.h"

#include "surfacenets.h"
#include "persistence.h"

const std::string path_model = "../tests/resources/Bear_2.ubs";

TEST(SurfaceNetsSuite, EstimateBounds) {
    const auto dao = persistence::MeshDao();
    const auto sdfn = dao.load_sdfn_from_file(path_model)[0];

    const auto result = surfacenets::estimate_bounding_box(sdfn);

    EXPECT_LE(result.row(0).sum(), 0);
    EXPECT_GE(0, result.row(0).sum());
}

TEST(SurfaceNetsSuite, E2E) {
    const auto dao = persistence::MeshDao();
    const auto sdfn = dao.load_sdfn_from_file(path_model)[0];

    const auto sut = surfacenets::SurfaceNetsMeshStrategy();

    const entities::QuadMesh result = sut.mesh(sdfn, 1000);

    EXPECT_EQ(0, result.n_vertices());
}
