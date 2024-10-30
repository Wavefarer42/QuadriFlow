#include <cstdlib>
#include <chrono>
#include <format>

#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include "OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh"
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

    field.initialize_parameterizer(args.face_count);

    spdlog::info("Elapsed: {:.3} seconds\n", watch);

    if (field.flag_preserve_boundary) {
        printf("Add boundary constrains...\n");
        Hierarchy &mRes = field.m_hierarchy;
        mRes.clearConstraints();
        for (uint32_t i = 0; i < 3 * mRes.mF.cols(); ++i) {
            if (mRes.mE2E[i] == -1) {
                uint32_t i0 = mRes.mF(i % 3, i / 3);
                uint32_t i1 = mRes.mF((i + 1) % 3, i / 3);
                Vector3d p0 = mRes.mV[0].col(i0), p1 = mRes.mV[0].col(i1);
                Vector3d edge = p1 - p0;
                if (edge.squaredNorm() > 0) {
                    edge.normalize();
                    mRes.mCO[0].col(i0) = p0;
                    mRes.mCO[0].col(i1) = p1;
                    mRes.mCQ[0].col(i0) = mRes.mCQ[0].col(i1) = edge;
                    mRes.mCQw[0][i0] = mRes.mCQw[0][i1] = mRes.mCOw[0][i0] = mRes.mCOw[0][i1] =
                            1.0;
                }
            }
        }
        mRes.propagateConstraints();
    }

    watch.reset();
    spdlog::info("Solving orientation field");

    Optimizer::optimize_orientations(field.m_hierarchy);
    field.find_orientation_singularities();

    spdlog::info("Elapsed: {:.3} seconds\n", watch);


    if (field.flag_adaptive_scale == 1) {
        watch.reset();
        spdlog::info("Analyzing mesh for adaptive scaling");

        field.estimate_slope();

        spdlog::info("Elapsed: {:.3} seconds\n", watch);
    }

    watch.reset();
    spdlog::info("Solving for field for adaptive scale");

    Optimizer::optimize_scale(field.m_hierarchy, field.rho, field.flag_adaptive_scale);
    field.flag_adaptive_scale = 1;

    spdlog::info("Elapsed: {:.3} seconds\n", watch);

    watch.reset();
    spdlog::info("Solving for position field");

    Optimizer::optimize_positions(field.m_hierarchy, field.flag_adaptive_scale);
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
