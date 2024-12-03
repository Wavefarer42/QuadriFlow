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


TEST(E2E, FullPipelineBoxComplex) {
    const auto dao = bootstrap::Container().mesh_dao();

    const auto dir_input = "../tests/resources/benchmark/7-box-sharpx-roundedy-rotated.ubs";
    const auto dir_output = "../tests/out/e2e";

    const auto faces = 10000;
    const auto path = fs::path(dir_input);;
    const auto path_base = std::format("{}/{}", dir_output, path.stem().string());

    if (fs::exists(path_base)) fs::remove_all(path_base);
    fs::create_directories(path_base);

    auto model = dao.load_model(path.string());
    auto sdfn = model[0];

    auto mesh = meshing::mesh_to_quadmesh(sdfn, model.bounding_box(0));
    dao.save_mesh(std::format("{}/{}.ply", path_base, "0-mesh"), mesh);

    mesh = smoothing::sdfn_projection(sdfn, mesh, 10);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "1-project"), mesh);

    mesh = smoothing::edge_snapping(sdfn, mesh, 10);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "2-intersect"), mesh);

    mesh = smoothing::laplacian_with_sdfn_projection(sdfn, mesh, 10);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "3-laplacian-project"), mesh);

    mesh = meshing::remesh_to_trimesh(mesh);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "4-trimesh"), mesh);

    mesh = meshing::remesh_to_quadmesh(sdfn, mesh, faces);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "5-remesh"), mesh);

    mesh = smoothing::sdfn_projection(sdfn, mesh);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "6-project"), mesh);

    mesh = smoothing::edge_snapping(sdfn, mesh);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "7-intersect"), mesh);

    // Notes on the pipeline
    // No final smooth as it breaks the edges.
}

TEST(E2E, Benchmark) {
    const auto service = bootstrap::Container().mesh_service();

    const auto dir_input = "../tests/resources/benchmark";
    const auto dir_output = "../tests/out/benchmark";

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
        const auto path_model = entry.path();
        const auto filename = path_model.filename().string();
        const auto ext = path_model.extension().string();

        for (int it_face: faces) {
            spdlog::info("--- Case {}/{} ---\n- filename:{}\n- faces: {}",
                         i_case, total_case, filename, it_face);
            const auto path_output = std::format("{}/{}", dir_output, filename);
            if (fs::exists(path_output)) {
                i_case++;
                continue;
            }

            try {
                service.to_isotropic_quadmesh(path_model, path_output, it_face);
            } catch (const std::exception &e) {
                spdlog::error("Case {}/{} failed: {}", i_case, total_case, e.what());
                failed++;
            }

            i_case++;
        }
    }

    spdlog::info("Finished benchmark with {} / {} failed cases", failed, total_case);
}
