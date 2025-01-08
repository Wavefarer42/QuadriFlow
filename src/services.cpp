#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

#include "services.h"
#include "meshing.h"
#include "smoothing.h"

namespace services {
    void MeshService::to_isotropic_quadmesh(
        const std::string &path_input,
        const std::string &path_output,
        const int face_count,
        const int sdfn_resolution
    ) const {
        spdlog::info(
            "Meshing unbound collection via isotropic quad meshing \n\tinput={}, output={}, faces={}, resolution={}",
            path_input, path_output, face_count, sdfn_resolution
        );

        std::filesystem::path _path_input = path_input;
        std::filesystem::path _path_output = path_output;

        assert(std::filesystem::exists(_path_input));
        assert(_path_input.extension() == ".ubs");
        assert(!std::filesystem::is_directory(_path_output) || std::filesystem::is_directory(_path_output));

        if (!std::filesystem::exists(_path_output)) {
            std::filesystem::create_directories(_path_output);
        }

        spdlog::stopwatch watch_total;

        entities::UnboundModel model = mesh_dao.load_model(path_input);
        for (int idx_model = 0; idx_model < model.size(); idx_model++) {
            try {
                spdlog::stopwatch watch;

                spdlog::info("Meshing model {}/{} via dual flow meshing", idx_model + 1, model.size());

                auto sdfn = model[idx_model];

                auto mesh = meshing::mesh_to_trimesh(sdfn, model.bounding_box(idx_model), sdfn_resolution);
                mesh = smoothing::fill_holes(mesh);
                mesh = smoothing::laplacian_with_sdfn_projection(sdfn, mesh, 20, 1);
                mesh = meshing::remesh_to_quadmesh(sdfn, mesh, face_count);
                mesh = smoothing::fill_holes(mesh);
                mesh = smoothing::sdfn_projection(sdfn, mesh);

                const auto path_filename = fmt::format(
                    "{}-{}-{}-{}.ply",
                    _path_input.stem().string(),
                    face_count,
                    sdfn_resolution,
                    idx_model
                );
                mesh_dao.save_mesh(_path_output / path_filename, mesh);

                spdlog::info("Finished meshing model {}/{} via dual flow meshing ({:.3}s)",
                             idx_model + 1, model.size(), watch);
            } catch (const std::exception &e) {
                spdlog::error("Error meshing model {}/{} via dual flow meshing: {}", idx_model + 1, model.size(),
                              e.what());
            }
        }

        spdlog::info("Finished meshing {} models via dual flow meshing faces={}, resolution={} ({:.3}s)",
                     model.size(), face_count, sdfn_resolution, watch_total);
    }
}
