#include "gtest/gtest.h"
#include <filesystem>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "sdfn.h"
#include "bootstrap.h"

namespace fs = std::filesystem;

TEST(E2E, FromSDFsSphere) {
    const auto service = bootstrap::Container().mesh_service();

    const Eigen::AlignedBox3f bounds(Eigen::Vector3f(-1.1, -1.1, -1.1), Eigen::Vector3f(1.1, 1.1, 1.1));
    auto mesh_quad = service.mesh_to_irregular_quadmesh(sdfn::sphere, bounds);
    const auto mesh_tri = service.remesh_to_trimesh(mesh_quad);
    const auto mesh_final = service.remesh_to_regular_quadmesh(mesh_tri);

    service.save_mesh("../tests/out/sphere-100-fff-10000-remesh.obj", mesh_final);
}

TEST(E2E, FromSDFsBox) {
    const auto service = bootstrap::Container().mesh_service();

    const Eigen::AlignedBox3f bounds(Eigen::Vector3f(-1.1, -1.1, -1.1), Eigen::Vector3f(1.1, 1.1, 1.1));
    auto mesh_quad = service.mesh_to_irregular_quadmesh(sdfn::box, bounds);
    const auto mesh_tri = service.remesh_to_trimesh(mesh_quad);
    const auto mesh_final = service.remesh_to_regular_quadmesh(mesh_tri);

    service.save_mesh("../tests/out/box-100-fff-10000-remesh.obj", mesh_final);
}

TEST(E2E, FromSDFsCylinder) {
    const auto service = bootstrap::Container().mesh_service();

    const Eigen::AlignedBox3f bounds(Eigen::Vector3f(-1.1, -1.1, -1.1), Eigen::Vector3f(1.1, 1.1, 1.1));
    auto mesh_quad = service.mesh_to_irregular_quadmesh(sdfn::cylinder, bounds);
    const auto mesh_tri = service.remesh_to_trimesh(mesh_quad);
    const auto mesh_final = service.remesh_to_regular_quadmesh(mesh_tri);

    service.save_mesh("../tests/out/cylinder-100-fff-10000-remesh.obj", mesh_final);
}

TEST(E2E, Benchmark) {
    const auto service = bootstrap::Container().mesh_service();

    const auto dir_input = "../tests/resources/benchmark/";
    const auto dir_output = "../tests/out/benchmark/";

    if (fs::exists(dir_output)) fs::remove_all(dir_output);
    fs::create_directories(dir_output);

    EXPECT_TRUE(fs::exists(dir_input));

    const std::vector faces = {100, 1000, 10000};
    const std::vector edges = {false, true};

    for (const auto &entry: fs::directory_iterator(dir_input)) {
        const auto path = entry.path();
        const auto filename = path.filename().string();
        const auto ext = path.extension().string();

        if (ext != ".ubs") continue;

        int total_case = faces.size() * edges.size();
        int i_case = 0;
        for (int it_face: faces) {
            for (bool it_edge: edges) {
                spdlog::info("--- Case {}/{}: {} ---", i_case, total_case, filename);

                auto model = service.load_unbound_model_from_file(path.string());
                auto mesh = service.mesh_to_irregular_quadmesh(model[0], model.bounding_box(0));
                mesh = service.remesh_to_trimesh(mesh);
                mesh = service.remesh_to_regular_quadmesh(
                    mesh, it_face, it_edge, false
                );

                const auto path_out = std::format("{}{}-f{}-e{}.ply",
                                                  dir_output, filename, it_face, it_edge ? "y" : "n");
                service.save_mesh(path_out, mesh);

                i_case++;
            }
        }
    }
}
