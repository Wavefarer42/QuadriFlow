#include "gtest/gtest.h"

#include "bootstrap.h"
#include "surfacenets.h"
#include "persistence.h"

const std::string path_model = "../tests/resources/Bear_2.ubs";

VectorXf sphere(const MatrixXf domain) {
    float radius = 1;
    const Vector3f origin = Vector3f::Zero();
    return (domain.rowwise() - origin.transpose()).rowwise().norm().array() - radius;
}

TEST(SurfaceNetsSuite, sphere) {
    MatrixXf domain(3, 3);
    domain << 0, 0, 0,
            1, 0, 0,
            0, 1, 0;
    const VectorXf result = sphere(domain);

    EXPECT_LE(result[0], 1);
    EXPECT_LE(result[1], 0);
    EXPECT_LE(result[2], 0);
}

TEST(SurfaceNetsSuite, EstimateBounds) {
    const auto sdfn = bootstrap::Container()
            .mesh_dao()
            .load_unbound_model(path_model)
            .sdfn_as_list()[0];

    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const AlignedBox3f result = sut.estimate_bounding_box(sphere, 20);

    EXPECT_EQ(result.min().x(), -1);
    EXPECT_EQ(result.min().y(), -1);
    EXPECT_EQ(result.min().z(), -1);
    EXPECT_EQ(result.max().x(), 1);
    EXPECT_EQ(result.max().y(), 1);
    EXPECT_EQ(result.max().z(), 1);
}

TEST(SurfaceNetsSuite, mesh) {
    const int resolution = 32;
    const AlignedBox3f bounds(Vector3f(-1.1, -1.1, -1.1), Vector3f(1.1, 1.1, 1.1));
    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(sphere, bounds, resolution);

    EXPECT_EQ(result.n_vertices(), 3992);
}

TEST(SurfaceNetsSuite, linear_indexing) {
    const auto resolution = 32;

    auto index = 0;
    MatrixXi indices(resolution * resolution * resolution, 3);
    for (int z = 0; z < resolution; ++z) {
        for (int y = 0; y < resolution; ++y) {
            for (int x = 0; x < resolution; ++x) {
                indices.row(index++) << x, y, z;
            }
        }
    }

    const auto linearize = [](Vector3i idx_nd) {
        return idx_nd.x() + idx_nd.y() * resolution + idx_nd.z() * resolution * resolution;
    };

    for (int i = 0; i < indices.rows(); ++i) {
        const int linear_index = linearize(indices.row(i));
        EXPECT_EQ(i, linear_index);
    }
}