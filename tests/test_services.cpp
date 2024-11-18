#include "gtest/gtest.h"
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "bootstrap.h"

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

    OpenMesh::IO::write_mesh(result, "../tests/out/box-gradient_smoothing.ply");
}

