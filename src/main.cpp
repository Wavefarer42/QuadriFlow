#include <filesystem>

#include <argparse/argparse.hpp>

#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

#include "bootstrap.h"
#include "services.h"

namespace fs = std::filesystem;
static const std::string version = "0.1.0";

struct CLIArgs {
    std::string path_in;
    std::string path_out;
    int face_count = -1;
    int preserve_edges = 0;
    int preserve_boundaries = 0;
    int use_adaptive_meshing = 0;
    int seed = 0;
    int resolution = 100;
    bool is_valid = true;
};


CLIArgs read_args(int argc, char **argv) {
    argparse::ArgumentParser program("meshbound");
    program.add_argument("--input")
            .help("Input unbound model file [.ubs]")
            .required();
    program.add_argument("--output")
            .help("Output directory.")
            .default_value(".");
    program.add_argument("--faces")
            .help("Target face count")
            .default_value(10000)
            .scan<'i', int>();
    program.add_argument("--resolution")
            .help("Resolution of the SDFn pre-meshing.")
            .default_value(100)
            .scan<'i', int>();


    int valid_args = 0;
    CLIArgs args;
    try {
        program.parse_args(argc, argv);
        args.path_in = program.get<std::string>("--input");
        args.path_out = program.get<std::string>("--output");
        args.face_count = program.get<int>("--faces");
        args.resolution = program.get<int>("--resolution");

        if (std::filesystem::exists(args.path_in) && args.path_in.ends_with(".ubs")) {
            valid_args++;
        } else {
            spdlog::error("Input file must exist and be .obj, .ply, or .ubs");
        }

        if (!fs::exists(args.path_out) || fs::is_directory(args.path_out)) {
            valid_args++;
        } else {
            spdlog::error("Output must be a directory or not exist");
        }

        if (args.face_count > 0) {
            valid_args++;
        } else {
            spdlog::error("Face count must be greater than 0");
        }

        if (args.resolution > 0) {
            valid_args++;
        } else {
            spdlog::error("Resolution must be greater than 0");
        }
    } catch (const std::runtime_error &err) {
        spdlog::error("Error parsing arguments: {}", err.what());
        spdlog::info("meshbound ({})\n{}", version, program.help().str());
    }

    if (valid_args == 4) {
        args.is_valid = true;
    } else {
        args.is_valid = false;
    }

    return args;
}


int main(int argc, char **argv) {
    const auto args = read_args(argc, argv);
    if (!args.is_valid) return 1;

    if (!fs::exists(args.path_out)) {
        fs::create_directories(args.path_out);
    }

    const spdlog::stopwatch watch_total;

    bootstrap::Container container = bootstrap::Container();
    const services::MeshService service = container.mesh_service();

    service.to_isotropic_quadmesh(
        args.path_in,
        args.path_out,
        args.resolution,
        args.face_count
    );

    spdlog::info("Finished generating mesh ({:.3}s)", watch_total);
    return 0;
}
