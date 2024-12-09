#include "gtest/gtest.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "Eigen/Core"

#include "sdfn.h"
#include "bootstrap.h"
#include "mathext.h"
#include "meshing.h"
#include "smoothing.h"

using namespace Eigen;


TEST(SmoothingTestSuite, SmoothingSDFnProjectionBox) {
    const auto dao = bootstrap::Container().mesh_dao();
    auto model = dao.load_model("../tests/resources/box-sharp-aligned.ubs");
    const auto sdfn = model[0];

    auto mesh_base = meshing::mesh_to_quadmesh(sdfn, model.bounding_box(0));
    auto mesh_smooth = mesh_base;
    const auto result = smoothing::sdfn_projection(sdfn, mesh_smooth, 10);

    OpenMesh::IO::write_mesh(mesh_base, "../tests/out/smoothing-sdfn-projection-box-base.ply");
    OpenMesh::IO::write_mesh(mesh_smooth, "../tests/out/smoothing-sdfn-projection-box-smooth.ply");
}

TEST(SmoothingTestSuite, LaplacianAngleFieldBox) {
    const AlignedBox3f bounds(Vector3f(-1.1, -1.1, -1.1), Vector3f(1.1, 1.1, 1.1));
    auto mesh = meshing::mesh_to_quadmesh(sdfn::box, bounds);
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

TEST(SmoothingTestSuite, LaplacianAngleFieldBoxComplex) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_model("../tests/resources/box-complex.ubs");

    auto mesh = meshing::mesh_to_quadmesh(model[0], model.bounding_box(0));
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

TEST(SmoothingTestSuite, SmoothingEdgeSnappingBox) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_model("../tests/resources/box-sharp-aligned.ubs");
    auto sdfn = model[0];

    auto mesh_base = meshing::mesh_to_quadmesh(sdfn, model.bounding_box(0));
    auto mesh_smooth = mesh_base;
    const auto result = smoothing::edge_snapping(sdfn, mesh_smooth, 10);

    OpenMesh::IO::write_mesh(mesh_base, "../tests/out/smoothing-edge-snapping-box-base.ply");
    OpenMesh::IO::write_mesh(mesh_smooth, "../tests/out/smoothing-edge-snapping-box-smooth.ply");
}

TEST(SmoothingTestSuite, SmoothingEdgeSnappingBoxSharpRoundedRotated) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_model("../tests/resources/benchmark-base/7-box-sharpx-roundedy-rotated.ubs");
    auto sdfn = model[0];

    auto mesh_base = meshing::mesh_to_quadmesh(sdfn, model.bounding_box(0));
    mesh_base = meshing::remesh_to_trimesh(mesh_base);
    auto mesh_smooth = mesh_base;
    const auto result = smoothing::edge_snapping(sdfn, mesh_smooth, 5);

    OpenMesh::IO::write_mesh(mesh_base, "../tests/out/smoothing-edge-snapping-box-base.ply");
    OpenMesh::IO::write_mesh(mesh_smooth, "../tests/out/smoothing-edge-snapping-box-smooth.ply");
}

TEST(SmoothingTestSuite, SmoothingLaplacianSDFnProjectionBoxSharpAligned) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_model("../tests/resources/box-sharp-aligned.ubs");
    auto sdfn = model[0];

    auto mesh_base = meshing::mesh_to_quadmesh(sdfn, model.bounding_box(0));
    auto mesh_smooth = mesh_base;
    const auto result = smoothing::laplacian_with_sdfn_projection(sdfn, mesh_smooth, 10);

    OpenMesh::IO::write_mesh(mesh_base, "../tests/out/smoothing-laplacian-projection-box-base.ply");
    OpenMesh::IO::write_mesh(mesh_smooth, "../tests/out/smoothing-laplacian-projection-box-smooth.ply");
}
