#include <filesystem>
#include <fstream>
#include <algorithm>
#include <Eigen/Core>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "gtest/gtest.h"

#include "bootstrap.h"
#include "smoothing.h"
#include "meshing.h"

namespace fs = std::filesystem;


TEST(E2E, FullPipelineStepwiseDebug) {
    const auto dao = bootstrap::Container().mesh_dao();

    const auto dir_input = "../tests/resources/benchmark/14-box-sphere-union-sharp.ubs";
    const auto dir_output = "../tests/out/e2e";

    const auto faces = 500;
    const auto faces_out = 100;
    const auto path = fs::path(dir_input);;
    const auto path_base = std::format("{}/{}", dir_output, path.stem().string());

    if (fs::exists(path_base)) fs::remove_all(path_base);
    fs::create_directories(path_base);

    auto model = dao.load_model(path.string());
    auto sdfn = model[0];

    auto mesh = meshing::mesh_to_trimesh(sdfn, model.bounding_box(0), 100);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "0-mesh"), mesh);

    mesh = smoothing::laplacian_with_sdfn_projection(sdfn, mesh, 10);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "1-laplacian-project"), mesh);

    mesh = smoothing::edge_snapping(sdfn, mesh, 10);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "2-project-edges"), mesh);

    mesh = meshing::remesh_to_quadmesh(sdfn, mesh, 20000);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "3-remesh-quad"), mesh);

    mesh = smoothing::sdfn_projection(sdfn, mesh);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "4-project"), mesh);

    mesh = smoothing::edge_snapping(sdfn, mesh);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "5-project-edges"), mesh);

    mesh = meshing::remesh_to_trimesh(mesh);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "6-remesh-tri"), mesh);

    mesh = meshing::remesh_to_quadmesh(sdfn, mesh, faces);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "7-remesh-quad"), mesh);

    mesh = smoothing::sdfn_projection(sdfn, mesh);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "8-project"), mesh);
}

TEST(E2E, Benchmark) {
    const auto service = bootstrap::Container().mesh_service();

    const auto dir_input = "../tests/resources/benchmark";
    const auto dir_output = "../tests/out/benchmark";

    if (!fs::exists(dir_output)) fs::create_directories(dir_output);

    EXPECT_TRUE(fs::exists(dir_input));

    const int sdfn_resolution = 100;
    const std::vector faces = {100, 1000, 10000};
    std::vector<fs::directory_entry> entries;
    for (const auto &entry: fs::directory_iterator(dir_input)) {
        if (entry.path().extension() == ".ubs" && !entry.path().string().ends_with(".skip.ubs")) {
            entries.emplace_back(entry);;
        }
    }
    std::ranges::sort(entries);

    int failed = 0;
    int i_case = 0;
    int total_case = faces.size() * entries.size();
    for (const auto &entry: entries) {
        const auto path_model = entry.path();

        for (int it_face: faces) {
            spdlog::info("--- Case {}/{} ---\n- filename:{}\n- faces: {}",
                         i_case, total_case, path_model.stem().string(), it_face);
            const auto path_output = std::format("{}/{}", dir_output, path_model.stem().string());
            const auto file_output = std::format("{}/{}-{}-{}-{}.ply",
                                                 path_output, path_model.stem().string(), it_face,
                                                 sdfn_resolution, 0);
            if (fs::exists(file_output)) {
                spdlog::info("Skipping case as output already exists");
                i_case++;
                continue;
            }

            try {
                service.to_isotropic_quadmesh(path_model, path_output, it_face, sdfn_resolution);
            } catch (const std::exception &e) {
                spdlog::error("Case {}/{} failed: {}", i_case, total_case, e.what());
                failed++;
            }

            i_case++;
        }
    }

    spdlog::info("Finished benchmark with {} / {} failed cases", failed, total_case);
}

TEST(E2E, PipelineDebug) {
    const auto service = bootstrap::Container().mesh_service();

    const fs::path model = "18-box-sphere-cut-roundedx-rotated.ubs";
    const fs::path path_input = "../tests/resources/benchmark" / model;
    const fs::path path_output = "../tests/out/benchmark" / path_input.stem();

    assert(fs::exists(path_input));
    if (fs::exists(path_output)) {
        fs::remove_all(path_output);
    }
    fs::create_directories(path_output);

    const std::vector param_resolution = {
        // 50,
        100,
        // 200
    };
    const std::vector param_face_count = {
        // 100,
        // 1000,
        10000
    };

    int failed = 0;
    int i_case = 0;
    int total_case = param_resolution.size() * param_face_count.size();

    for (int it_resolution: param_resolution) {
        for (int it_face: param_face_count) {
            spdlog::info(
                "--- Case {}/{} ---\n- filename:{}\n- faces: {} \n- resolution: {}",
                i_case,
                total_case,
                path_input.stem().string(),
                it_face,
                it_resolution
            );

            try {
                service.to_isotropic_quadmesh(path_input, path_output, it_face, it_resolution);
            } catch (const std::exception &e) {
                spdlog::error("Case {}/{} failed: {}", i_case, total_case, e.what());
                failed++;
            }

            i_case++;
        }
    }


    spdlog::info("Finished benchmark with {} / {} failed cases", failed, total_case);
}
