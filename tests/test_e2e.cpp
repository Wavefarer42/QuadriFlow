#include "gtest/gtest.h"
#include <filesystem>
#include <fstream>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "sdfn.h"
#include "bootstrap.h"
#include "smoothing.h"

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

TEST(E2E, FullPipeline) {
    const auto service = bootstrap::Container().mesh_service();

    const auto dir_input = "../tests/resources/benchmark/";
    const auto dir_output = "../tests/out/e2e";


    EXPECT_TRUE(fs::exists(dir_input));

    auto i = 0;
    const auto faces = 10000;
    for (const auto &entry: fs::directory_iterator(dir_input)) {
        if (i > 0) break;
        const auto path = entry.path();
        const auto path_base = std::format("{}/{}", dir_output, path.stem().string());

        if (fs::exists(path_base)) fs::remove_all(path_base);
        fs::create_directories(path_base);

        if (path.extension().string() != ".ubs") continue;

        auto model = service.load_unbound_model_from_file(path.string());
        auto sdfn = model[0];

        auto mesh = service.mesh_to_irregular_quadmesh(sdfn, model.bounding_box(0));
        service.save_mesh(std::format("{}/{}.ply", path_base, "0-mesh"), mesh);

        mesh = smoothing::sdfn_projection(sdfn, mesh, 10);
        service.save_mesh(std::format("{}/{}.ply", path_base, "1-project"), mesh);

        mesh = smoothing::edge_snapping(sdfn, mesh, 10);
        service.save_mesh(std::format("{}/{}.ply", path_base, "2-intersect"), mesh);

        mesh = smoothing::laplacian_with_sdfn_projection(sdfn, mesh, 10);
        service.save_mesh(std::format("{}/{}.ply", path_base, "3-laplacian-project"), mesh);

        mesh = service.remesh_to_trimesh(mesh);
        service.save_mesh(std::format("{}/{}.ply", path_base, "4-trimesh"), mesh);

        mesh = service.remesh_to_regular_quadmesh(mesh, faces, true, false);
        service.save_mesh(std::format("{}/{}.ply", path_base, "5-remesh"), mesh);

        mesh = smoothing::sdfn_projection(sdfn, mesh);
        service.save_mesh(std::format("{}/{}.ply", path_base, "6-project"), mesh);

        mesh = smoothing::edge_snapping(sdfn, mesh);
        service.save_mesh(std::format("{}/{}.ply", path_base, "7-intersect"), mesh);

        // Notes on the pipeline
        // No final smooth as it breaks the edges.
        i++;
    }
}

TEST(E2E, FullPipelineBoxComplex) {
    const auto service = bootstrap::Container().mesh_service();

    const auto dir_input = "../tests/resources/benchmark/7-box-sharpx-roundedy-rotated.ubs";
    const auto dir_output = "../tests/out/e2e";

    const auto faces = 10000;
    const auto path = fs::path(dir_input);;
    const auto path_base = std::format("{}/{}", dir_output, path.stem().string());

    if (fs::exists(path_base)) fs::remove_all(path_base);
    fs::create_directories(path_base);

    auto model = service.load_unbound_model_from_file(path.string());
    auto sdfn = model[0];

    auto mesh = service.mesh_to_irregular_quadmesh(sdfn, model.bounding_box(0));
    service.save_mesh(std::format("{}/{}.ply", path_base, "0-mesh"), mesh);

    mesh = service.remesh_to_trimesh(mesh);
    service.save_mesh(std::format("{}/{}.ply", path_base, "1-trimesh"), mesh);

    mesh = smoothing::laplacian_with_sdfn_projection(sdfn, mesh, 10);
    service.save_mesh(std::format("{}/{}.ply", path_base, "2-laplacian-project"), mesh);

    mesh = service.remesh_to_regular_quadmesh(mesh, faces, true, false);
    service.save_mesh(std::format("{}/{}.ply", path_base, "3-remesh"), mesh);

    mesh = smoothing::sdfn_projection(sdfn, mesh);
    service.save_mesh(std::format("{}/{}.ply", path_base, "4-project"), mesh);

    // Need to fix the degeneration of the mesh
    mesh = smoothing::edge_snapping(sdfn, mesh);
    service.save_mesh(std::format("{}/{}.ply", path_base, "5-intersect"), mesh);

    mesh = service.remesh_to_trimesh(mesh);
    service.save_mesh(std::format("{}/{}.ply", path_base, "6-trimesh"), mesh);

    mesh = service.remesh_to_regular_quadmesh(mesh, faces, true, false);
    service.save_mesh(std::format("{}/{}.ply", path_base, "7-remesh"), mesh);


    // Notes on the pipeline
    // No final smooth as it breaks the edges.
}

TEST(E2E, Benchmark) {
    const auto service = bootstrap::Container().mesh_service();

    const auto dir_input = "../tests/resources/benchmark/";
    const auto dir_output = "../tests/out/benchmark/";

    if (!fs::exists(dir_output)) fs::create_directories(dir_output);

    EXPECT_TRUE(fs::exists(dir_input));

    const std::vector faces = {100, 1000, 10000};
    std::vector<fs::directory_entry> entries;
    for (const auto &entry: fs::directory_iterator(dir_input)) {
        if (entry.path().extension() == ".ubs" && !entry.path().string().ends_with(".skip.ubs")) {
            entries.emplace_back(entry);;
        }
    }
    std::sort(entries.begin(), entries.end());

    int failed = 0;
    int i_case = 0;
    int total_case = faces.size() * entries.size();
    for (const auto &entry: entries) {
        const auto path = entry.path();
        const auto filename = path.filename().string();
        const auto ext = path.extension().string();

        for (int it_face: faces) {
            spdlog::info("--- Case {}/{} ---\n- filename:{}\n- faces: {}",
                         i_case, total_case, filename, it_face);
            const auto path_out = std::format("{}{}-f{}.ply",
                                              dir_output, filename, it_face);
            if (fs::exists(path_out)) {
                i_case++;
                continue;
            }

            try {
                auto model = service.load_unbound_model_from_file(path.string());
                auto sdfn = model[0];

                auto mesh = service.mesh_to_irregular_quadmesh(sdfn, model.bounding_box(0));
                mesh = service.remesh_to_trimesh(mesh);
                mesh = smoothing::laplacian_with_sdfn_projection(sdfn, mesh, 10);
                mesh = service.remesh_to_regular_quadmesh(mesh, it_face, true);
                mesh = smoothing::sdfn_projection(sdfn, mesh);

                service.save_mesh(path_out, mesh);
            } catch (const std::exception &e) {
                spdlog::error("Case {}/{} failed: {}", i_case, total_case, e.what());
                failed++;
            }

            i_case++;
        }
    }

    spdlog::info("Finished benchmark with {} / {} failed cases", failed, total_case);
}
