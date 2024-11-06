#include <Eigen/Dense>

#include "surfacenets.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

using namespace Eigen;

namespace surfacenets {

    MatrixXf cartesian_product(
            const VectorXf &X,
            const VectorXf &Y,
            const VectorXf &Z
    ) {
        spdlog::debug("Computing cartesian product");

        MatrixXf domain(X.size() * Y.size() * Z.size(), 3);
        int i = 0;
        for (float x: X) {
            for (float y: Y) {
                for (float z: Z) {
                    domain.row(i++) << x, y, z;
                }
            }
        }

        return domain;
    }

    MatrixXf create_sample_grid(
            const AlignedBox3f &bounds,
            float resolution
    ) {
        spdlog::debug("Creating sample grid");

        const auto X = VectorXf::LinSpaced(static_cast<int>(resolution), bounds.min().x(), bounds.max().x());
        const auto Y = VectorXf::LinSpaced(static_cast<int>(resolution), bounds.min().y(), bounds.max().y());
        const auto Z = VectorXf::LinSpaced(static_cast<int>(resolution), bounds.min().z(), bounds.max().z());

        return cartesian_product(X, Y, Z);
    }

    AlignedBox3f estimate_bounding_box(
            entities::SDFn sdfn,
            float initial_resolution,
            float threshold = 1e-3
    ) {
        spdlog::debug("Estimating bounding box");
        spdlog::stopwatch watch;

        // Start with a reasonably sized initial bounding box
        AlignedBox3f bounds(Vector3f(-10.0f, -10.0f, -10.0f),
                            Vector3f(10.0f, 10.0f, 10.0f));
        AlignedBox3f bounds_prev;

        // Set the maximum number of iterations and adaptive resolution settings
        const int max_rounds = 32;
        float resolution = initial_resolution;
        const float resolution_decay = 0.5f;  // Reduce resolution for finer sampling

        int round = 0;
        while (round < max_rounds) {
            // Create the sample grid within the current bounding box
            const auto domain = create_sample_grid(bounds, resolution);
            const auto sdf = sdfn(domain);

            // Initialize new bounds to be updated based on the SDF evaluation
            AlignedBox3f new_bounds = bounds;

            for (auto i = 0; i < domain.rows(); ++i) {
                const auto distance = std::abs(sdf(i));
                if (distance > threshold) {
                    const VectorXf &point = domain.row(i);
                    new_bounds.min() = new_bounds.min().cwiseMin(point);
                    new_bounds.max() = new_bounds.max().cwiseMax(point);
                }
            }

            // Check for convergence by comparing the volume change
            float volume_change = (new_bounds.max() - new_bounds.min()).norm() -
                                  (bounds.max() - bounds.min()).norm();
            if (std::abs(volume_change) < 1e-3) {
                break;
            }

            // Update the bounds and reduce resolution for finer iteration
            bounds = new_bounds;
            resolution *= resolution_decay;
            round++;
        }

        spdlog::debug("Estimated bounding box: min = ({}, {}, {}), max = ({}, {}, {}), rounds = {}, elapsed = {}",
                      bounds.min().x(), bounds.min().y(), bounds.min().z(),
                      bounds.max().x(), bounds.max().y(), bounds.max().z(),
                      round, watch);
        return bounds;
    }

    AlignedBox3f estimate_bounding_box(
            entities::SDFn sdfn,
            int grid_resolution
    ) {
        spdlog::debug("Estimating bounding box");
        spdlog::stopwatch watch;

        float step_size = 0.5f;

        VectorXf lower_bound = VectorXf::Constant(3, std::numeric_limits<float>::max());
        VectorXf upper_bound = VectorXf::Constant(3, std::numeric_limits<float>::lowest());

        const int total_points = (2 * grid_resolution + 1) * (2 * grid_resolution + 1) * (2 * grid_resolution + 1);

        MatrixXf domain(total_points, 3);
        int index = 0;
        for (int x = -grid_resolution; x <= grid_resolution; ++x) {
            for (int y = -grid_resolution; y <= grid_resolution; ++y) {
                for (int z = -grid_resolution; z <= grid_resolution; ++z) {
                    domain.row(index++) << x * step_size, y * step_size, z * step_size;
                }
            }
        }

        const VectorXf distances = sdfn(domain);

        for (int i = 0; i < total_points; ++i) {
            if (distances[i] <= 0.0f) {
                lower_bound = lower_bound.array().min(domain.row(i).transpose().array());
                upper_bound = upper_bound.array().max(domain.row(i).transpose().array());
            }
        }

        AlignedBox3f bounds(lower_bound, upper_bound);
        spdlog::debug("Estimated bounding box: min = ({}, {}, {}), max = ({}, {}, {}), elapsed = {}",
                      bounds.min().x(), bounds.min().y(), bounds.min().z(),
                      bounds.max().x(), bounds.max().y(), bounds.max().z(), watch);
        return bounds;
    }

    entities::QuadMesh estimate_surface_vertices(
            const entities::SDFn sdfn,
            const MatrixXf &domain,
            const MatrixXf &sdf
    ) {
        spdlog::debug("Estimating surface vertices");
        entities::QuadMesh mesh;

        for (int i = 0; i < sdf.rows(); ++i) {
        }

        return mesh;
    }

    entities::QuadMesh create_surface_faces(
            const entities::SDFn sdfn,
            const MatrixXf &domain,
            const entities::QuadMesh &mesh
    ) {
        spdlog::debug("Creating surface faces");

        return mesh;
    }

    entities::QuadMesh SurfaceNetsMeshStrategy::mesh(const entities::SDFn sdfn, const int resolution) const {
        spdlog::debug("Meshing SDFn with resolution {}", resolution);

        const auto bounds = estimate_bounding_box(sdfn, resolution);
        const auto domain = create_sample_grid(bounds, resolution);
        const MatrixXf sdf = sdfn(domain);
        entities::QuadMesh mesh = estimate_surface_vertices(sdfn, domain, sdf);
        mesh = create_surface_faces(sdfn, domain, mesh);

        return mesh;
    }

}