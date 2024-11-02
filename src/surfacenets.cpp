#include <Eigen/Dense>

#include "surfacenets.h"
#include "spdlog/spdlog.h"

using namespace Eigen;

namespace surfacenets {

    MatrixXf cartesian_product(const ArrayXf X, const ArrayXf Y, const ArrayXf Z) {
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

    MatrixXf estimate_bounding_box(const entities::SDFn sdfn) {
        spdlog::debug("Estimating bounding box");

        const float step = 16;
        float x0 = -1e9, y0 = -1e9, z0 = -1e9;
        float x1 = 1e9, y1 = 1e9, z1 = 1e9;

        float threshold_prev = 1e9;
        for (int i = 0; i < 32; ++i) {
            const auto X = ArrayXf::LinSpaced(step, x0, x1);
            const auto Y = ArrayXf::LinSpaced(step, y0, y1);
            const auto Z = ArrayXf::LinSpaced(step, z0, z1);
            const Vector3f distance(X[1] - X[0], Y[1] - Y[0], Z[1] - Z[0]);
            const auto threshold = distance.norm() / 2;
            if (std::abs(threshold - threshold_prev) < 1e-4) {
                break;
            }

            threshold_prev = threshold;
            const auto domain = cartesian_product(X, Y, Z);
            const VectorXf sdf = sdfn(domain).col(0);

            int lenX = X.size(), lenY = Y.size(), lenZ = Z.size();
            MatrixXd volumeMatrix(lenX, lenY * lenZ);
            for (int idx = 0; idx < sdf.size(); ++idx) {
                int x_idx = idx / (lenY * lenZ);
                int yz_idx = idx % (lenY * lenZ);
                volumeMatrix(x_idx, yz_idx) = sdf(idx);
            }

            std::vector<Vector3i> withinThreshold;
            for (int idx = 0; idx < sdf.size(); ++idx) {
                if (std::abs(sdf(idx)) <= threshold) {
                    int z_idx = idx % lenZ;
                    int y_idx = (idx / lenZ) % lenY;
                    int x_idx = idx / (lenY * lenZ);
                    withinThreshold.emplace_back(x_idx, y_idx, z_idx);
                }
            }

            Vector3i where_max = withinThreshold[0];
            Vector3i where_min = withinThreshold[0];
            for (const auto &vec: withinThreshold) {
                where_max = where_max.cwiseMax(vec);
                where_min = where_min.cwiseMin(vec);
            }

            x1 = x0 + where_max(0) * distance(0) + distance(0) / 2;
            y1 = y0 + where_max(1) * distance(1) + distance(1) / 2;
            z1 = z0 + where_max(2) * distance(2) + distance(2) / 2;

            x0 = x0 + where_min(0) * distance(0) - distance(0) / 2;
            y0 = y0 + where_min(1) * distance(1) - distance(1) / 2;
            z0 = z0 + where_min(2) * distance(2) - distance(2) / 2;

        }

        MatrixXf result(2, 3);
        result << x0, y0, z0, x1, y1, z1;
        return result;
    }

    MatrixXf create_sample_grid(
            const MatrixXf &bounds,
            const int resolution
    ) {
        spdlog::debug("Creating sample grid");

        const auto X = ArrayXf::LinSpaced(resolution, bounds(0, 0), bounds(1, 0));
        const auto Y = ArrayXf::LinSpaced(resolution, bounds(0, 1), bounds(1, 1));
        const auto Z = ArrayXf::LinSpaced(resolution, bounds(0, 2), bounds(1, 2));

        return cartesian_product(X, Y, Z);
    }

    entities::QuadMesh estimate_surface_vertices(
            const entities::SDFn sdfn,
            const MatrixXf &domain,
            const MatrixXf &sdf
    ) {
        spdlog::debug("Estimating surface vertices");
        entities::QuadMesh mesh;

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

        const auto bounds = estimate_bounding_box(sdfn);
        const auto domain = create_sample_grid(bounds, resolution);
        const MatrixXf sdf = sdfn(domain);
        entities::QuadMesh mesh = estimate_surface_vertices(sdfn, domain, sdf);
        mesh = create_surface_faces(sdfn, domain, mesh);

        return mesh;
    }

}