#include <filesystem>

#include <argparse/argparse.hpp>

#include "adapters.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

#include "bootstrap.h"
#include "services.h"
#include "field-math.h"

using namespace services;
const std::string version = "0.1.0";

entities::CLIArgs read_args(int argc, char **argv) {
    argparse::ArgumentParser program("meshbound");
    program.add_argument("--input")
            .help("Input mesh file [obj, ply, ubs]")
            .required();
    program.add_argument("--output")
            .help("Output mesh file [obj, ply]")
            .required();
    program.add_argument("--edges")
            .help("Detect and preserve sharp edges")
            .default_value(false)
            .implicit_value(true);
    program.add_argument("--boundary")
            .help("Preserve shape boundaries")
            .default_value(false)
            .implicit_value(true);
    program.add_argument("--adaptive")
            .help("Adaptive scaling")
            .default_value(false)
            .implicit_value(true);
    program.add_argument("--seed")
            .help("Seed of the run")
            .default_value(14)
            .scan<'i', int>();
    program.add_argument("--faces")
            .help("Target face count")
            .default_value(10000)
            .scan<'i', int>();
    program.add_argument("--resolution")
            .help("Resolution of the SDFn pre-meshing.")
            .default_value(100)
            .scan<'i', int>();

    spdlog::info("meshbound ({})\n{}", version, program.help().str());

    int valid_args = 0;
    entities::CLIArgs args;
    try {
        program.parse_args(argc, argv);
        args.path_in = program.get<std::string>("--input");
        args.path_out = program.get<std::string>("--output");
        args.face_count = program.get<int>("--faces");
        args.use_adaptive_meshing = program.get<bool>("--adaptive");
        args.preserve_boundaries = program.get<bool>("--boundary");
        args.preserve_edges = program.get<bool>("--edges");
        args.seed = program.get<int>("--seed");
        args.resolution = program.get<int>("--resolution");

        if (std::filesystem::exists(args.path_in)
            && (args.path_in.ends_with(".obj")
                || args.path_in.ends_with(".ply")
                || args.path_in.ends_with(".ubs"))) {
            valid_args++;
        } else {
            spdlog::error("Input file must exist and be .obj, .ply, or .ubs");
        }

        if (args.path_out.ends_with(".obj") || args.path_out.ends_with(".ply")) {
            valid_args++;
        } else {
            spdlog::error("Output file must be .obj or .ply");
        }

        if (args.face_count > 0) {
            valid_args++;
        } else {
            spdlog::error("Face count must be greater than 0");
        }

        if (args.seed > 0) {
            valid_args++;
        } else {
            spdlog::error("Seed must be greater than 0");
        }

        if (args.resolution > 0) {
            valid_args++;
        } else {
            spdlog::error("Resolution must be greater than 0");
        }
    } catch (const std::runtime_error &err) {
        spdlog::error("Error parsing arguments: {}", err.what());
    }

    if (valid_args == 5) {
        args.is_valid = true;
    } else {
        args.is_valid = false;
    }

    return args;
}


int main(int argc, char **argv) {
    spdlog::stopwatch watch_total;

    const auto args = read_args(argc, argv);
    if (!args.is_valid) return 1;

    bootstrap::Container container = bootstrap::Container();
    const MeshService service = container.mesh_service();

    entities::Mesh mesh;
    std::optional<std::reference_wrapper<entities::SDFn> > model_opt = std::nullopt;
    if (args.path_in.ends_with(".ubs")) {
        auto model = service.load_unbound_model_from_file(args.path_in);
        if (model.size() > 0) {
            mesh = service.mesh_to_irregular_quadmesh(model[0], model.bounding_box(0), args.resolution);
            model_opt = model[0];
        } else {
            spdlog::warn("The given Unbound collection does not contain any models.");
            return 2;
        }
    } else {
        mesh = service.load_mesh(args.path_in);
    }

    service.remesh_to_trimesh(mesh);

    services::Parametrizer field = service.remesh_to_regular_quadmesh(
        mesh,
        args.face_count,
        args.preserve_edges,
        args.use_adaptive_meshing,
        model_opt
    );

    // service.save_mesh(args.path_out, field);
    service.save_mesh(args.path_out, adapters::from_parametrizer_to_quad_mesh(field));

    spdlog::info("Finished generating mesh ({:.3}s)", watch_total);
    return 0;
}
