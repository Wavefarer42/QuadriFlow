#include "gtest/gtest.h"

#include "sdfn.h"
#include "bootstrap.h"
#include <Eigen/Core>
#include <Eigen/Dense>

TEST(E2E, FromSDFsSphere) {
    const auto service = bootstrap::Container().mesh_service();

    const int resolution = 100;
    const Eigen::AlignedBox3f bounds(Eigen::Vector3f(-1.1, -1.1, -1.1), Eigen::Vector3f(1.1, 1.1, 1.1));
    auto mesh_quad = service.mesh(sdfn::sphere, resolution, bounds);
    const auto mesh_tri = service.to_trimesh(mesh_quad);
    const auto mesh_final = service.remesh(mesh_tri);

    service.save_mesh("../tests/out/sphere-100-fff-10000-remesh.obj", mesh_final);
}

TEST(E2E, FromSDFsBox) {
    const auto service = bootstrap::Container().mesh_service();

    const int resolution = 100;
    const Eigen::AlignedBox3f bounds(Eigen::Vector3f(-1.1, -1.1, -1.1), Eigen::Vector3f(1.1, 1.1, 1.1));
    auto mesh_quad = service.mesh(sdfn::box, resolution, bounds);
    const auto mesh_tri = service.to_trimesh(mesh_quad);
    const auto mesh_final = service.remesh(mesh_tri);

    service.save_mesh("../tests/out/box-100-fff-10000-remesh.obj", mesh_final);
}

TEST(E2E, FromSDFsCylinder) {
    const auto service = bootstrap::Container().mesh_service();

    const int resolution = 100;
    const Eigen::AlignedBox3f bounds(Eigen::Vector3f(-1.1, -1.1, -1.1), Eigen::Vector3f(1.1, 1.1, 1.1));
    auto mesh_quad = service.mesh(sdfn::cylinder, resolution, bounds);
    const auto mesh_tri = service.to_trimesh(mesh_quad);
    const auto mesh_final = service.remesh(mesh_tri);

    service.save_mesh("../tests/out/cylinder-100-fff-10000-remesh.obj", mesh_final);
}