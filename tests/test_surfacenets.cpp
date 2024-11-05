#include "gtest/gtest.h"

#include "bootstrap.h"
#include "surfacenets.h"
#include "persistence.h"

const std::string path_model = "../tests/resources/Bear_2.ubs";

TEST(SurfaceNetsSuite, EstimateBounds) {
    const auto sdfn = bootstrap::Container()
            .mesh_dao()
            .load_unbound_model(path_model)
            .sdfn_as_list()[0];

    const AlignedBox3f result = surfacenets::estimate_bounding_box(sdfn, 100, 1e-3);

    EXPECT_LE(result.min().x(), 0);
    EXPECT_LE(result.min().y(), 0);
    EXPECT_LE(result.min().z(), 0);
    EXPECT_GE(result.max().x(), 0);
    EXPECT_GE(result.max().y(), 0);
    EXPECT_GE(result.max().z(), 0);
}
