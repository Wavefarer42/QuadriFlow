#include <filesystem>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <regex>
#include <string>
#include <Eigen/Core>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "gtest/gtest.h"

#include "bootstrap.h"
#include "smoothing.h"
#include "meshing.h"

namespace fs = std::filesystem;

TEST(E2E, FullPipelineStepwiseDebug) {
    const auto dao = bootstrap::Container().mesh_dao();

    const auto dir_input = "../tests/resources/benchmark/A-01.only-19.ubs";
    // const auto dir_input = "../tests/resources/benchmark/A-02.only-14.ubs";
    // const auto dir_input = "../tests/resources/benchmark/A-06.only.ubs";
    // const auto dir_input = "../tests/resources/sausage.ubs";
    const auto dir_output = "../tests/out/benchmark";

    const auto idx_model = 19;
    const auto faces = 10000;
    const auto path = fs::path(dir_input);
    const auto path_base = std::format("{}/{}", dir_output, path.stem().string());

    if (fs::exists(path_base)) fs::remove_all(path_base);
    fs::create_directories(path_base);

    auto model = dao.load_model(path.string());
    auto sdfn = model[idx_model];

    auto mesh = meshing::mesh_to_trimesh(sdfn, model.bounding_box(idx_model), 200, "marching-cubes-33");
    dao.save_mesh(std::format("{}/{}.ply", path_base, "0-mesh"), mesh);

    mesh = smoothing::laplacian_with_sdfn_projection(sdfn, mesh, 30);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "1-laplacian-project"), mesh);

    mesh = smoothing::edge_snapping(sdfn, mesh, 3);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "2-edges"), mesh);

    mesh = meshing::remesh_to_quadmesh(sdfn, mesh, faces);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "3-remesh-quad"), mesh);

    mesh = smoothing::sdfn_projection(sdfn, mesh, 5);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "4-project"), mesh);

    mesh = smoothing::laplacian_with_sdfn_projection(sdfn, mesh, 20);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "5-project"), mesh);

    mesh = smoothing::edge_snapping(sdfn, mesh, 3);
    dao.save_mesh(std::format("{}/{}.ply", path_base, "6-edges"), mesh);
}


TEST(E2E, Benchmark) {
    const std::vector param_resolution = {
        // 50,
        // 100,
        200
    };
    const std::vector param_face_count = {
        // 100,
        // 1000,
        // 5000,
        10000
    };

    const fs::path path_input = "../tests/resources/benchmark";
    const fs::path path_output = "../tests/out/benchmark";

    assert(fs::exists(path_input));
    assert(param_resolution.size() > 0);
    assert(param_face_count.size() > 0);

    const std::regex re_onlies(R"(\.only\.|\.only-\d+)");
    const std::regex re_only_number(R"(\.only-(\d+))");
    std::vector<fs::directory_entry> onlies;
    std::vector<fs::directory_entry> paths_collections;
    for (const auto &entry: fs::directory_iterator(path_input)) {
        if (entry.path().extension() == ".ubs" && !entry.path().string().ends_with(".skip.ubs")) {
            paths_collections.emplace_back(entry);
            if (entry.path().filename().string().find(".only") != std::string::npos) {
                onlies.emplace_back(entry);
            }
        }
    }
    if (!onlies.empty()) paths_collections = onlies;
    std::ranges::sort(paths_collections);

    const auto dao = bootstrap::Container().mesh_dao();

    fs::create_directories(path_output);

    int failed = 0;
    for (auto &it_col: paths_collections) {
        const auto filename = it_col.path().filename().string();
        auto model = dao.load_model(it_col.path().string());

        int model_only = -1;
        std::smatch match;
        if (std::regex_search(filename, match, re_only_number)) {
            if (match.size() > 1) {
                model_only = std::stoi(match[1]);
            }
        }

        for (int idx_model = 0; idx_model < model.size(); ++idx_model) {
            if (model_only >= 0 && model_only != idx_model) continue;

            for (int it_resolution: param_resolution) {
                for (int it_face: param_face_count) {
                    spdlog::info(
                        "-------- Collection: {}\nmodel: {}/{}\nresolution: {}\nfaces: {}",
                        it_col.path().stem().string(),
                        idx_model,
                        model.size(),
                        it_resolution,
                        it_face
                    );

                    const std::string path_base = std::format(
                        "{}/{}-{}-{}-{}",
                        path_output.string(), it_col.path().stem().string(),
                        idx_model, it_resolution, it_face
                    );
                    if (fs::exists(path_base)) continue;

                    fs::create_directories(path_base);

                    try {
                        auto sdfn = model[idx_model];

                        auto mesh = meshing::mesh_to_trimesh(sdfn, model.bounding_box(idx_model), it_resolution,
                                                             "marching-cubes-33");
                        dao.save_mesh(std::format("{}/{}.ply", path_base, "0-mesh"), mesh);

                        mesh = smoothing::laplacian_with_sdfn_projection(sdfn, mesh, 30);
                        dao.save_mesh(std::format("{}/{}.ply", path_base, "1-laplacian-project"), mesh);

                        mesh = smoothing::edge_snapping(sdfn, mesh, 3);
                        dao.save_mesh(std::format("{}/{}.ply", path_base, "2-edges"), mesh);

                        mesh = meshing::remesh_to_quadmesh(sdfn, mesh, it_face);
                        dao.save_mesh(std::format("{}/{}.ply", path_base, "3-remesh-quad"), mesh);

                        mesh = smoothing::sdfn_projection(sdfn, mesh, 5);
                        dao.save_mesh(std::format("{}/{}.ply", path_base, "4-project"), mesh);
                    } catch (const std::exception &e) {
                        spdlog::error(
                            "-------- Collection: {}\nmodel: {}\nresolution: {}\nfaces: {}\nerror: {}",
                            it_col.path().stem().string(),
                            idx_model,
                            it_resolution,
                            it_face,
                            e.what()
                        );
                        std::filesystem::path log_file = std::format("{}/{}", path_base, "failed.txt");
                        std::ofstream log_stream(log_file, std::ios::app);
                        log_stream << "-------- Collection: " << it_col.path().stem().string() << "\n"
                                << "model: " << idx_model << "\n"
                                << "resolution: " << it_resolution << "\n"
                                << "faces: " << it_face << "\n"
                                << "error: " << e.what() << "\n"
                                << "----------------------\n";
                        failed++;
                    }catch (...) {
                        spdlog::error(
                            "-------- Collection: {}\nmodel: {}\nresolution: {}\nfaces: {}",
                            it_col.path().stem().string(),
                            idx_model,
                            it_resolution,
                            it_face
                        );
                        failed++;
                    }
                }
            }
        }
    }

    spdlog::info("Finished benchmark with {} failed cases", failed);
}
