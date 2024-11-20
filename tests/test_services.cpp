#include "gtest/gtest.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "Eigen/Core"

#include "sdfn.h"
#include "bootstrap.h"

using namespace Eigen;

TEST(MeshServiceSuite, to_trimesh) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const auto service = bootstrap::Container().mesh_service();
    auto mesh = service.load_mesh("../tests/resources/sphere.ply");

    const auto result = service.to_trimesh(mesh);

    OpenMesh::IO::write_mesh(result, "../tests/out/sphere-to_trimesh.ply");

    EXPECT_TRUE(service.is_trimesh(mesh));
}

TEST(MeshServiceSuite, GradientSmoothingBox) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const auto service = bootstrap::Container().mesh_service();
    auto mesh = service.load_mesh("../tests/resources/box.ply");

    const auto result = service.to_trimesh(mesh);

    OpenMesh::IO::write_mesh(result, "../tests/out/box-smoothing_surface_snapping.ply");
}

TEST(MeshServiceSuite, LaplacianAngleField) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const auto service = bootstrap::Container().mesh_service();

    const int resolution = 100;
    const AlignedBox3f bounds(Vector3f(-1.1, -1.1, -1.1), Vector3f(1.1, 1.1, 1.1));
    auto mesh = service.mesh(sdfn::box, resolution, bounds);
    const Vector3f result = service.create_laplacian_angle_field(sdfn::box, mesh);

    const float threshold = result.mean() + (2 * (result.array() - result.mean()).square().sum() / (result.size() - 1));

    for (auto it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
        mesh.delete_face(*it, false);
    }

    for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
        if (result[it->idx()] < threshold) {
            mesh.delete_vertex(*it, false);
        }
    }

    mesh.garbage_collection();
    OpenMesh::IO::write_mesh(mesh, "../tests/out/box-laplacian_angle_field.ply");
}