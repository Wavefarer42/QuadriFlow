#include <cstdlib>
#include <format>
#include <filesystem>

#include <argparse/argparse.hpp>
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

#include "adapters.h"
#include "bootstrap.h"
#include "services.h"
#include "optimizer.h"
#include "field-math.h"

using namespace services;

Parametrizer field;

entities::CLIArgs read_args(int argc, char **argv) {

    argparse::ArgumentParser program("meshbound");
    program.add_argument("--input")
            .help("Input mesh file [obj, ply, ubs]")
            .required();
    program.add_argument("--output")
            .help("Output mesh file [obj, ply]")
            .required();
    program.add_argument("--sharp")
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

    spdlog::info("meshbound \n{}", program.help().str());

    int valid_args = 0;
    entities::CLIArgs args;
    try {
        program.parse_args(argc, argv);
        args.path_in = program.get<std::string>("--input");
        args.path_out = program.get<std::string>("--output");
        args.face_count = program.get<int>("--faces");
        args.adaptive = program.get<bool>("--adaptive");
        args.boundaries = program.get<bool>("--boundary");
        args.edges = program.get<bool>("--sharp");
        args.seed = program.get<int>("--seed");

        if (std::filesystem::exists(args.path_in)
            && (args.path_in.ends_with(".obj")
                || args.path_in.ends_with(".ply"))) {
            valid_args++;
        } else {
            spdlog::error("Input file must exist and be .obj or .ply");
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
    } catch (const std::runtime_error &err) {
        spdlog::error("Error parsing arguments: {}", err.what());
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


    bootstrap::Container container = bootstrap::Container();
    const MeshService service = container.mesh_service();
    const auto mesh = service.load_trimesh_from_file(args.path_in);
    adapters::initialize_parameterizer(field, mesh);

    spdlog::stopwatch watch, watch_total;
    spdlog::info("Initializing parameters");

    field.initialize_parameterizer(
            args.boundaries,
            args.edges,
            args.face_count,
            args.adaptive
    );

    spdlog::info("Elapsed: {:.3} seconds\n", watch);

    if (args.boundaries) {
        service.set_boundary_constraints(field.m_hierarchy);
    }

    watch.reset();
    spdlog::info("Solving orientation field");

    Optimizer::optimize_orientations(field.m_hierarchy);
    service.find_orientation_singularities(field.m_hierarchy);

    spdlog::info("Elapsed: {:.3} seconds\n", watch);


    if (args.adaptive) {
        watch.reset();
        spdlog::info("Analyzing mesh for adaptive scaling");

        const auto [faces_slope, faces_orientation] = service.estimate_slope(
                field.m_hierarchy,
                field.m_triangle_space,
                field.m_faces_normals
        );
        field.m_faces_slope = faces_slope;
        field.m_faces_orientation = faces_orientation;

        spdlog::info("Elapsed: {:.3} seconds\n", watch);
    }

    watch.reset();
    spdlog::info("Solving field for adaptive scale");

    Optimizer::optimize_scale(field.m_hierarchy, field.m_rho, args.adaptive);

    spdlog::info("Elapsed: {:.3} seconds\n", watch);

    watch.reset();
    spdlog::info("Solving for position field");

    Optimizer::optimize_positions(field.m_hierarchy);
    const auto [singularity_position, singularity_rank, singularity_index] = service.find_position_singularities(
            field.m_hierarchy,
            true
    );
    field.m_singularity_position = singularity_position;
    field.m_singularity_rank = singularity_rank;
    field.m_singularity_index = singularity_index;

    spdlog::info("Elapsed: {:.3} seconds\n", watch);

    watch.reset();
    spdlog::info("Solving for integer constraints");

    field.compute_index_map(field.m_hierarchy);

    spdlog::info("Elapsed: {:.3} seconds\n", watch);

    service.save_quadmesh_to_file(args.path_out, field);

    spdlog::info("Total elapsed: {:.3} seconds\n", watch_total);
    return 0;
}
