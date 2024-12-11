#include "gtest/gtest.h"
#include "spdlog/spdlog.h"
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "sdfn.h"
#include "bootstrap.h"
#include "meshing.h"

using namespace Eigen;

TEST(MeshingTestSuite, SphereSDF) {
    MatrixXf domain(3, 3);
    domain << 0, 0, 0,
            1, 0, 0,
            0, 1, 0;
    const VectorXf result = sdfn::sphere(domain);

    EXPECT_LE(result[0], 1);
    EXPECT_LE(result[1], 0);
    EXPECT_LE(result[2], 0);
}

TEST(MeshingTestSuite, MeshingSphere) {
    spdlog::set_level(spdlog::level::debug);

    const AlignedBox3f bounds(Vector3f(-1.1, -1.1, -1.1), Vector3f(1.1, 1.1, 1.1));
    const auto result = meshing::mesh_to_quadmesh(sdfn::sphere, bounds, 100, "surface-nets");

    OpenMesh::IO::write_mesh(result, "../tests/out/sphere.ply");

    EXPECT_EQ(result.n_vertices(), 39008);
}

TEST(MeshingTestSuite, MeshingBox) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const AlignedBox3f bounds(Vector3f(-2, -2, -2), Vector3f(2, 2, 2));
    const auto result = meshing::mesh_to_quadmesh(sdfn::box, bounds, 100, "surface-nets");

    OpenMesh::IO::write_mesh(result, "../tests/out/box.ply");

    EXPECT_EQ(result.n_vertices(), 49688);
}

TEST(MeshingTestSuite, MeshingCylinder) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const AlignedBox3f bounds(Vector3f(-1.1, -1.1, -1.1), Vector3f(1.1, 1.1, 1.1));
    const auto result = meshing::mesh_to_quadmesh(sdfn::cylinder, bounds, 100, "surface-nets");

    OpenMesh::IO::write_mesh(result, "../tests/out/cylinder.ply");
}

TEST(MeshingTestSuite, MeshingRotatedBox) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const AlignedBox3f bounds(Vector3f(-5, -5, -5), Vector3f(5, 5, 5));
    const auto sdfn = sdfn::rotate(sdfn::box, Vector3f(0, 1, 0), 0.5);
    const auto result = meshing::mesh_to_quadmesh(sdfn, bounds, 100, "surface-nets");

    OpenMesh::IO::write_mesh(result, "../tests/out/box-rotated.ply");
}

TEST(MeshingTestSuite, to_trimesh) {
    const auto dao = bootstrap::Container().mesh_dao();
    auto mesh = dao.load_mesh("../tests/resources/sphere.ply");

    const auto result = meshing::remesh_to_trimesh(mesh);

    OpenMesh::IO::write_mesh(result, "../tests/out/sphere-to_trimesh.ply");

    EXPECT_TRUE(meshing::is_trimesh(mesh));
}

TEST(MeshingTestSuite, LinearIndexing) {
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

TEST(MeshingTestSuite, MeshingUnboundBear) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_model("../tests/resources/Bear_2.ubs");

    const auto result = meshing::mesh_to_quadmesh(model[0], model.bounding_box(0), 100, "surface-nets");

    OpenMesh::IO::write_mesh(result, "../tests/out/Bear_2.ply");

    EXPECT_EQ(result.n_vertices(), 25346);
}

TEST(MeshingTestSuite, MeshingUnboundBox) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_model("../tests/resources/box.ubs");

    const auto result = meshing::mesh_to_quadmesh(model[0], model.bounding_box(0), 100, "surface-nets");

    OpenMesh::IO::write_mesh(result, "../tests/out/box-unbound.ply");

    EXPECT_EQ(result.n_vertices(), 58568);
}

TEST(MeshingTestSuite, MeshingUnboundBoxComplex) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_model("../tests/resources/box-complex.ubs");

    const auto result = meshing::mesh_to_quadmesh(model[0], model.bounding_box(0), 100, "surface-nets");

    OpenMesh::IO::write_mesh(result, "../tests/out/box-complex.ply");

    EXPECT_EQ(result.n_vertices(), 34766);
}

TEST(MeshingTestSuite, MeshingUnboundBoxSharpAligned) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_model("../tests/resources/box-sharp-aligned.ubs");

    const auto result = meshing::mesh_to_quadmesh(model[0], model.bounding_box(0), 100, "surface-nets");

    OpenMesh::IO::write_mesh(result, "../tests/out/box-sharp-aligned.ply");

    EXPECT_EQ(result.n_vertices(), 58018);
}

TEST(MeshingTestSuite, MeshingUnboundBoxSharpRoundAligned) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_model("../tests/resources/box-sharp-round-aligned.ubs");

    const auto result = meshing::mesh_to_quadmesh(model[0], model.bounding_box(0), 100, "surface-nets");

    OpenMesh::IO::write_mesh(result, "../tests/out/box-sharp-round-aligned.ply");

    EXPECT_EQ(result.n_vertices(), 56770);
}

TEST(MeshingTestSuite, MeshingUnboundBoxSharpRoundRotated) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_model("../tests/resources/box-sharp-round-rotated.ubs");

    const auto result = meshing::mesh_to_quadmesh(model[0], model.bounding_box(0), 100, "surface-nets");

    OpenMesh::IO::write_mesh(result, "../tests/out/box-sharp-round-rotated.ply");

    EXPECT_EQ(result.n_vertices(), 0);
}

TEST(MeshingTestSuite, MeshingBoxSphereCut) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_model("../tests/resources/benchmark/08-case.skip.ubs");

    const auto result = meshing::mesh_to_quadmesh(model[0], model.bounding_box(0), 100, "surface-nets");

    OpenMesh::IO::write_mesh(result, "../tests/out/08-case.ply");

    EXPECT_EQ(result.n_vertices(), 0);
}

TEST(MeshingTestSuite, MeshingUnboundBoxSharpRoundRotatedDelaunay) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_model("../tests/resources/box-sharp-round-rotated.ubs");

    const auto result = meshing::mesh_to_trimesh(model[0], model.bounding_box(0), 100, "delaunay");

    OpenMesh::IO::write_mesh(result, "../tests/out/box-sharp-round-rotated.ply");

    EXPECT_EQ(result.n_vertices(), 0);
}
