#include "gtest/gtest.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "Eigen/Core"

#include "bootstrap.h"
#include "mathext.h"

using namespace Eigen;

TEST(MeshingSuite, RemeshNoEdgeConstraints) {
    const auto service = bootstrap::Container().mesh_service();
    auto model = service.load_unbound_model_from_file("../tests/resources/box-complex.ubs");
    auto mesh = service.mesh_to_irregular_quadmesh(model[0], model.bounding_box(0));
    mesh = service.remesh_to_trimesh(mesh);

    entities::Mesh result = service.remesh_to_regular_quadmesh(
        mesh, 10000, true, false
    );

    service.save_mesh("../tests/out/box-complex-remesh.ply", result);
}

TEST(MeshingSuite, RemeshEdgeConstraints) {
    const auto service = bootstrap::Container().mesh_service();
    auto model = service.load_unbound_model_from_file("../tests/resources/box-complex.ubs");
    auto mesh = service.mesh_to_irregular_quadmesh(model[0], model.bounding_box(0));
    mesh = service.remesh_to_trimesh(mesh);

    entities::Mesh result = service.remesh_to_regular_quadmesh(
        mesh, 10000, true, false, model[0]
    );

    service.save_mesh("../tests/out/box-complex-remesh-edges.ply", result);
}


TEST(MeshingSuite, RemeshSharpEdgesBoxComplex) {
    const auto service = bootstrap::Container().mesh_service();
    auto model = service.load_unbound_model_from_file("../tests/resources/box-complex.ubs");
    auto mesh = service.mesh_to_irregular_quadmesh(model[0], model.bounding_box(0));
    mesh = service.remesh_to_trimesh(mesh);

    entities::Mesh mesh_base = service.remesh_to_regular_quadmesh(
        mesh, 10000, false, false
    );

    entities::Mesh mesh_edges = service.remesh_to_regular_quadmesh(
        mesh, 10000, true, false, model[0]
    );

    service.save_mesh("../tests/out/box-complex-remesh-base.ply", mesh_base);
    service.save_mesh("../tests/out/box-complex-remesh-edges.ply", mesh_edges);
}
