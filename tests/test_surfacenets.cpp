#include "gtest/gtest.h"
#include "spdlog/spdlog.h"
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "sdfn.h"
#include "bootstrap.h"
#include "surfacenets.h"
#include "persistence.h"

const std::string path_model = "../tests/resources/Bear_2.ubs";


TEST(SurfaceNetsSuite, SphereSDF) {
    MatrixXf domain(3, 3);
    domain << 0, 0, 0,
            1, 0, 0,
            0, 1, 0;
    const VectorXf result = sdfn::sphere(domain);

    EXPECT_LE(result[0], 1);
    EXPECT_LE(result[1], 0);
    EXPECT_LE(result[2], 0);
}

TEST(SurfaceNetsSuite, EstimateBounds) {
    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const AlignedBox3f result = sut.estimate_bounding_box(sdfn::sphere, 20);

    EXPECT_EQ(result.min().x(), -1);
    EXPECT_EQ(result.min().y(), -1);
    EXPECT_EQ(result.min().z(), -1);
    EXPECT_EQ(result.max().x(), 1);
    EXPECT_EQ(result.max().y(), 1);
    EXPECT_EQ(result.max().z(), 1);
}

TEST(SurfaceNetsSuite, Meshing) {
    spdlog::set_level(spdlog::level::debug);

    const int resolution = 20;
    const AlignedBox3f bounds(Vector3f(-1.1, -1.1, -1.1), Vector3f(1.1, 1.1, 1.1));
    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(sdfn::sphere, resolution, bounds);

    OpenMesh::IO::write_mesh(result, "../tests/out/sphere.ply");

    EXPECT_EQ(result.n_vertices(), 39008);
}

TEST(SurfaceNetsSuite, MeshingBox) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const int resolution = 30;
    const AlignedBox3f bounds(Vector3f(-1.1, -1.1, -1.1), Vector3f(1.1, 1.1, 1.1));
    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(sdfn::box, resolution, bounds);

    OpenMesh::IO::write_mesh(result, "../tests/out/box.ply");

    EXPECT_EQ(result.n_vertices(), 39008);
}

TEST(SurfaceNetsSuite, MeshingCylinder) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const int resolution = 100;
    const AlignedBox3f bounds(Vector3f(-1.1, -1.1, -1.1), Vector3f(1.1, 1.1, 1.1));
    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(sdfn::cylinder, resolution, bounds);

    OpenMesh::IO::write_mesh(result, "../tests/out/cylinder.ply");

    EXPECT_EQ(result.n_vertices(), 39008);
}

TEST(SurfaceNetsSuite, MeshingRotatedBox) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const int resolution = 100;
    const AlignedBox3f bounds(Vector3f(-5, -5, -5), Vector3f(5, 5, 5));
    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const auto sdfn = sdfn::rotate(sdfn::box, Vector3f(1, 0, 0), 0.5);
    const entities::Mesh result = sut.mesh(sdfn, resolution, bounds);

    OpenMesh::IO::write_mesh(result, "../tests/out/box-rotated.ply");

    EXPECT_EQ(result.n_vertices(), 39008);
}

TEST(SurfaceNetsSuite, LinearIndexing) {
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

TEST(SurfaceNetsSuite, MeshingUnboundModel) {
    spdlog::set_level(spdlog::level::debug);

    const auto sdfn = bootstrap::Container()
            .mesh_dao()
            .load_unbound_model(path_model)
            .sdfn_as_list()[1];

    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(sdfn);

    OpenMesh::IO::write_mesh(result, "../tests/out/bear.ply");

    EXPECT_EQ(result.n_vertices(), 1928);
}
