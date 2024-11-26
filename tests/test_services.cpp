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

TEST(MeshServiceSuite, GradientSmoothingBoxComplex) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const auto service = bootstrap::Container().mesh_service();
    const auto sdfn = service.load_unbound_model_from_file("../tests/resources/box-complex.ubs")[0];
    const int resolution = 100;
    const AlignedBox3f bounds(Vector3f(-200, -200, -200), Vector3f(200, 200, 200));

    auto mesh = service.mesh(sdfn, resolution, bounds);

    const auto result = service.smoothing_surface_snapping(sdfn, mesh, 10);
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

    mesh.request_face_status();
    for (auto it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
        mesh.delete_face(*it, false);
    }

    mesh.request_vertex_status();
    for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
        if (field[it->idx()] < 30) {
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

    mesh.request_face_status();
    for (auto it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
        mesh.delete_face(*it, false);
    }

    mesh.request_vertex_status();
    for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
        if (field[it->idx()] < 30) {
            mesh.delete_vertex(*it, false);
        }
    }

    mesh.garbage_collection();
    OpenMesh::IO::write_mesh(mesh, "../tests/out/box-complex-laplacian_angle_field.ply");
}

TEST(MeshServiceSuite, SmoothingEdgeSnappingUnboundBox) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const auto service = bootstrap::Container().mesh_service();
    const auto sdfn = service.load_unbound_model_from_file("../tests/resources/box.ubs")[0];
    const int resolution = 100;
    const AlignedBox3f bounds(Vector3f(-200, -200, -200), Vector3f(200, 200, 200));

    auto mesh = service.mesh(sdfn, resolution, bounds);
    const auto result = service.smoothing_edge_snapping(sdfn, mesh, 3);

    OpenMesh::IO::write_mesh(result, "../tests/out/box-smoothing_edge_snapping.ply");
}

TEST(MeshServiceSuite, SmoothingEdgeSnappingUnboundBoxComplex) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const auto service = bootstrap::Container().mesh_service();
    const auto sdfn = service.load_unbound_model_from_file("../tests/resources/box-complex.ubs")[0];
    const int resolution = 100;
    const AlignedBox3f bounds(Vector3f(-200, -200, -200), Vector3f(200, 200, 200));

    auto mesh = service.mesh(sdfn, resolution, bounds);
    const auto result = service.smoothing_edge_snapping(sdfn, mesh, 1);

    OpenMesh::IO::write_mesh(result, "../tests/out/box-complex-smoothing_edge_snapping.ply");
}

TEST(MeshServiceSuite, SmoothingLaplacianSDFnProjectionBoxComplex) {
#ifdef DEV_DEBUG
    spdlog::set_level(spdlog::level::debug);
#endif

    const auto service = bootstrap::Container().mesh_service();
    const auto sdfn = service.load_unbound_model_from_file("../tests/resources/box-complex.ubs")[0];
    const int resolution = 100;
    const AlignedBox3f bounds(Vector3f(-200, -200, -200), Vector3f(200, 200, 200));

    auto mesh = service.mesh(sdfn, resolution, bounds);
    const auto result = service.smoothing_laplacian_sdf_projection(sdfn, mesh, 10);

    OpenMesh::IO::write_mesh(result, "../tests/out/smoothing-laplacian-sdfn-projection-box-complex.ply");
}
