#include "gtest/gtest.h"

#include "sdfn.h"
#include "bootstrap.h"
#include <Eigen/Core>
#include <Eigen/Dense>

TEST(E2E, FromSDFs) {
    const auto service = bootstrap::Container().mesh_service();

    const int resolution = 100;
    const Eigen::AlignedBox3f bounds(Eigen::Vector3f(-1.1, -1.1, -1.1), Eigen::Vector3f(1.1, 1.1, 1.1));
    const auto mesh = service.mesh(sdfn::sphere, bounds, resolution);
    const auto remesh = service.remesh(mesh);

    service.save_mesh("../tests/out/sphere-100-fff-10000-remesh.ply", remesh);

// FIXME mesh is quad mesh but a trimesh i expected
    EXPECT_EQ(1000, mesh.n_vertices());
}