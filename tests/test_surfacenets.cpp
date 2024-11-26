#include "gtest/gtest.h"
#include "spdlog/spdlog.h"
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "sdfn.h"
#include "bootstrap.h"
#include "surfacenets.h"


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

TEST(SurfaceNetsSuite, MeshingSphere) {
    spdlog::set_level(spdlog::level::debug);

    const AlignedBox3f bounds(Vector3f(-1.1, -1.1, -1.1), Vector3f(1.1, 1.1, 1.1));
    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(sdfn::sphere, bounds);

    OpenMesh::IO::write_mesh(result, "../tests/out/sphere.ply");

    EXPECT_EQ(result.n_vertices(), 1568);
}

TEST(SurfaceNetsSuite, MeshingBox) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const AlignedBox3f bounds(Vector3f(-1.1, -1.1, -1.1), Vector3f(1.1, 1.1, 1.1));
    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(sdfn::box, bounds);

    OpenMesh::IO::write_mesh(result, "../tests/out/box.ply");

    EXPECT_EQ(result.n_vertices(), 49688);
}

TEST(SurfaceNetsSuite, MeshingCylinder) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const AlignedBox3f bounds(Vector3f(-1.1, -1.1, -1.1), Vector3f(1.1, 1.1, 1.1));
    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(sdfn::cylinder, bounds);

    OpenMesh::IO::write_mesh(result, "../tests/out/cylinder.ply");
}

TEST(SurfaceNetsSuite, MeshingRotatedBox) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const AlignedBox3f bounds(Vector3f(-5, -5, -5), Vector3f(5, 5, 5));
    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const auto sdfn = sdfn::rotate(sdfn::box, Vector3f(0, 1, 0), 0.5);
    const entities::Mesh result = sut.mesh(sdfn, bounds);

    OpenMesh::IO::write_mesh(result, "../tests/out/box-rotated.ply");
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

TEST(SurfaceNetsSuite, MeshingUnboundBear) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_unbound_model("../tests/resources/Bear_2.ubs");

    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(model[0], model.bounding_box(0), 100);

    OpenMesh::IO::write_mesh(result, "../tests/out/Bear_2.ply");

    EXPECT_EQ(result.n_vertices(), 5677);
}

TEST(SurfaceNetsSuite, MeshingUnboundBox) {
    const auto sdfn = bootstrap::Container()
            .mesh_dao()
            .load_unbound_model("../tests/resources/box.ubs")[0];

    const AlignedBox3f bounds(Vector3f(-51, -51, -51), Vector3f(51, 51, 51));
    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(sdfn, bounds);

    OpenMesh::IO::write_mesh(result, "../tests/out/box-unbound.ply");

    EXPECT_EQ(result.n_vertices(), 58568);
}

TEST(SurfaceNetsSuite, MeshingUnboundBoxComplex) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_unbound_model("../tests/resources/box-complex.ubs");

    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(model[0], model.bounding_box(0));

    OpenMesh::IO::write_mesh(result, "../tests/out/box-complex.ply");

    EXPECT_EQ(result.n_vertices(), 6578);
}

TEST(SurfaceNetsSuite, MeshingUnboundBoxSharpAligned) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_unbound_model("../tests/resources/box-sharp-aligned.ubs");

    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(model[0], model.bounding_box(0));

    OpenMesh::IO::write_mesh(result, "../tests/out/box-sharp-aligned.ply");

    EXPECT_EQ(result.n_vertices(), 34766);
}

TEST(SurfaceNetsSuite, MeshingUnboundBoxSharpRoundAligned) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_unbound_model("../tests/resources/box-sharp-round-aligned.ubs");

    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(model[0], model.bounding_box(0));

    OpenMesh::IO::write_mesh(result, "../tests/out/box-sharp-round-aligned.ply");

    EXPECT_EQ(result.n_vertices(), 56770);
}

TEST(SurfaceNetsSuite, MeshingUnboundBoxSharpRoundRotated) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_unbound_model("../tests/resources/box-sharp-round-rotated.ubs");

    const auto sut = surfacenets::SurfaceNetsMeshStrategy();
    const entities::Mesh result = sut.mesh(model[0], model.bounding_box(0));

    OpenMesh::IO::write_mesh(result, "../tests/out/box-sharp-round-rotated.ply");

    EXPECT_EQ(result.n_vertices(), 5677);
}
