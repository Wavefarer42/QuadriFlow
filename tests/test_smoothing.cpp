#include "gtest/gtest.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "Eigen/Core"

#include "sdfn.h"
#include "bootstrap.h"
#include "mathext.h"
#include "smoothing.h"

using namespace Eigen;

TEST(MeshServiceSuite, to_trimesh) {
    const auto service = bootstrap::Container().mesh_service();
    auto mesh = service.load_mesh("../tests/resources/sphere.ply");

    const auto result = service.remesh_to_trimesh(mesh);

    OpenMesh::IO::write_mesh(result, "../tests/out/sphere-to_trimesh.ply");

    EXPECT_TRUE(service.is_trimesh(mesh));
}

TEST(MeshServiceSuite, GradientSmoothingBoxComplex) {
    const auto service = bootstrap::Container().mesh_service();
    auto model = service.load_unbound_model_from_file("../tests/resources/box-complex.ubs");

    auto mesh = service.mesh_to_irregular_quadmesh(model[0], model.bounding_box(0));

    const auto result = smoothing::sdfn_projection(model[0], mesh, 10);
}

TEST(MeshServiceSuite, LaplacianAngleFieldBox) {
    const auto service = bootstrap::Container().mesh_service();

    const AlignedBox3f bounds(Vector3f(-1.1, -1.1, -1.1), Vector3f(1.1, 1.1, 1.1));
    auto mesh = service.mesh_to_irregular_quadmesh(sdfn::box, bounds);
    const VectorXf field = smoothing::laplacian_angular_field(sdfn::box, mesh);

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
    const auto service = bootstrap::Container().mesh_service();
    auto model = service.load_unbound_model_from_file("../tests/resources/box-complex.ubs");

    auto mesh = service.mesh_to_irregular_quadmesh(model[0], model.bounding_box(0));
    const VectorXf field = smoothing::laplacian_angular_field(model[0], mesh);

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
    const auto service = bootstrap::Container().mesh_service();
    auto model = service.load_unbound_model_from_file("../tests/resources/box.ubs");

    auto mesh = service.mesh_to_irregular_quadmesh(model[0], model.bounding_box(0));
    const auto result = smoothing::edge_snapping(model[0], mesh, 3);

    OpenMesh::IO::write_mesh(result, "../tests/out/box-smoothing_edge_snapping.ply");
}

TEST(MeshServiceSuite, SmoothingEdgeSnappingUnboundBoxComplex) {
    const auto service = bootstrap::Container().mesh_service();
    auto model = service.load_unbound_model_from_file("../tests/resources/box-complex.ubs");

    auto mesh = service.mesh_to_irregular_quadmesh(model[0], model.bounding_box(0));
    const auto result = smoothing::edge_snapping(model[0], mesh, 1);

    OpenMesh::IO::write_mesh(result, "../tests/out/box-complex-smoothing_edge_snapping.ply");
}

TEST(MeshServiceSuite, SmoothingLaplacianSDFnProjectionBoxComplex) {
    const auto service = bootstrap::Container().mesh_service();
    auto model = service.load_unbound_model_from_file("../tests/resources/box-complex.ubs");

    auto mesh = service.mesh_to_irregular_quadmesh(model[0], model.bounding_box(0));
    const auto result = smoothing::laplacian_with_sdfn_projection(model[0], mesh, 10);

    OpenMesh::IO::write_mesh(result, "../tests/out/smoothing-laplacian-sdfn-projection-box-complex.ply");
}
