#include "gtest/gtest.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "Eigen/Core"

#include "sdfn.h"
#include "bootstrap.h"
#include "mathext.h"

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

TEST(MeshServiceSuite, LaplacianAngleFieldBox) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const auto service = bootstrap::Container().mesh_service();

    const int resolution = 100;
    const AlignedBox3f bounds(Vector3f(-1.1, -1.1, -1.1), Vector3f(1.1, 1.1, 1.1));
    auto mesh = service.mesh(sdfn::box, resolution, bounds);
    const VectorXf field = service.create_laplacian_angle_field(sdfn::box, mesh);

    std::cout << "Field: " << mathext::count_unique(field) << std::endl;

    const float threshold = field.mean();

    mesh.request_face_status();
    for (auto it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
        mesh.delete_face(*it, false);
    }

    mesh.request_vertex_status();
    for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
        if (field[it->idx()] < threshold) {
            mesh.delete_vertex(*it, false);
        }
    }

    mesh.garbage_collection();
    OpenMesh::IO::write_mesh(mesh, "../tests/out/box-laplacian_angle_field.ply");
}

TEST(MeshServiceSuite, LaplacianAngleFieldUnboundBoxComplex) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const auto service = bootstrap::Container().mesh_service();
    const auto sdfn = service.load_unbound_model_from_file("../tests/resources/box-complex.ubs")[0];
    const int resolution = 100;
    const AlignedBox3f bounds(Vector3f(-200, -200, -200), Vector3f(200, 200, 200));

    auto mesh = service.mesh(sdfn, resolution, bounds);
    const VectorXf field = service.create_laplacian_angle_field(sdfn, mesh);

    std::cout << "Field: " << mathext::count_unique(field) << std::endl;

    const float threshold = field.mean();

    mesh.request_face_status();
    for (auto it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
        mesh.delete_face(*it, false);
    }

    mesh.request_vertex_status();
    for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
        if (field[it->idx()] < threshold) {
            mesh.delete_vertex(*it, false);
        }
    }

    mesh.garbage_collection();
    OpenMesh::IO::write_mesh(mesh, "../tests/out/box-complex-laplacian_angle_field.ply");
}
