#include <Eigen/Dense>

#include "surfacenets.h"
#include "spdlog/spdlog.h"

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
            float resolution,
            float threshold = 1e-3
    ) {
        spdlog::debug("Estimating bounding box");

        AlignedBox3f bounds(Vector3f(-1e18, -1e18, -1e18),
                            Vector3f(1e18, 1e18, 1e18));
        AlignedBox3f bounds_prev(Vector3f(0, 0, 0),
                                 Vector3f(0, 0, 0));

        const int max_rounds = 32;
        int round = 0;
        while (round < max_rounds &&
               (bounds.max() != bounds_prev.max() || bounds.min() != bounds_prev.min())) {
            const auto domain = create_sample_grid(bounds, resolution);
            const auto sdf = sdfn(domain);

            for (auto i = 0; i < domain.rows(); ++i) {
                if (std::abs(sdf(i)) < threshold) {
                    const auto &point = domain.row(i);
                    if (point.x() <= bounds.min().x() + resolution) bounds.min().x() -= resolution;
                    if (point.x() >= bounds.max().x() - resolution) bounds.max().x() += resolution;
                    if (point.y() <= bounds.min().y() + resolution) bounds.min().y() -= resolution;
                    if (point.y() >= bounds.max().y() - resolution) bounds.max().y() += resolution;
                    if (point.z() <= bounds.min().z() + resolution) bounds.min().z() -= resolution;
                    if (point.z() >= bounds.max().z() - resolution) bounds.max().z() += resolution;
                }
            }

            bounds_prev = bounds;
            round++;
        }

        spdlog::debug("Estimated bounding box: min = ({}, {}, {}), max = ({}, {}, {}), rounds = {}",
                      bounds.min().x(), bounds.min().y(), bounds.min().z(),
                      bounds.max().x(), bounds.max().y(), bounds.max().z(),
                      round);
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