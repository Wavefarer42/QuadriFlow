#include <cstdlib>
#include <format>

#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include "adapters.h"
#include "bootstrap.h"
#include "services.h"
#include "optimizer.h"
#include "field-math.h"

using namespace services;

Parametrizer field;

unsigned long long inline GetCurrentTime64() {
    using namespace std::chrono;
    return duration_cast<milliseconds>(steady_clock::now().time_since_epoch()).count();
}

entities::CLIArgs read_args(int argc, char **argv) {

    entities::CLIArgs args;
    for (int i = 0; i < argc; ++i) {
        if (strcmp(argv[i], "--input") == 0) {
            args.path_in = argv[i + 1];
        } else if (strcmp(argv[i], "--output") == 0) {
            args.path_out = argv[i + 1];
        } else if (strcmp(argv[i], "--sharp") == 0) {
            args.edges = 1;
        } else if (strcmp(argv[i], "--boundary") == 0) {
            args.boundaries = 1;
        } else if (strcmp(argv[i], "--adaptive") == 0) {
            args.adaptive = 1;
        } else if (strcmp(argv[i], "--seed") == 0) {
            field.m_hierarchy.rng_seed = std::stoi(argv[i + 1]);
        } else if (strcmp(argv[i], "--faces") == 0) {
            args.face_count = std::stoi(argv[i + 1]);
        }
    }

    if (args.path_in.empty()) {
        spdlog::error("Expected input mesh file [obj, ply] with argument --input");
        exit(1);
    }

    if (args.path_out.empty()) {
        spdlog::error("Expected output mesh file [obj, ply] with argument --output");
        exit(1);
    }

    return args;
}


int main(int argc, char **argv) {
    spdlog::info("Meshbound --input <path>.[obj,ply] --output <path>.[obj,ply]");

    const auto args = read_args(argc, argv);

    bootstrap::Container container = bootstrap::Container();
    MeshService service = container.mesh_service();
    const auto mesh = service.load_trimesh_from_file(args.path_in);
    adapters::initialize_parameterizer(field, mesh);

    spdlog::stopwatch watch, watch_total;
    spdlog::info("Initializing parameters");

    field.initialize_parameterizer(args.face_count, args.adaptive);

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
    field.find_position_singularities();

    spdlog::info("Elapsed: {:.3} seconds\n", watch);

    watch.reset();
    spdlog::info("Solving for integer constraints");

    field.compute_index_map();

    spdlog::info("Elapsed: {:.3} seconds\n", watch);

    field.save_to_obj(args.path_out.c_str());

    spdlog::info("Total elapsed: {:.3} seconds\n", watch_total);
    return 0;
}
