#pragma once

#include "entities.h"

using namespace Eigen;

namespace surfacenets {
    using NdToFlatIndexer = std::function<int(const Vector3i &)>;

    class SurfaceNetsMeshStrategy {
    public:
        [[nodiscard]]
        entities::Mesh mesh(
                const entities::SDFn &sdfn,
                const AlignedBox3f &bounds,
                int resolution = 100
        ) const;

    private:
        const MatrixXi CUBE_CORNERS = (MatrixXi(8, 3) << 0, 0, 0,
                1, 0, 0,
                0, 1, 0,
                1, 1, 0,
                0, 0, 1,
                1, 0, 1,
                0, 1, 1,
                1, 1, 1).finished();

        const Matrix<int, 12, 2> CUBE_EDGES = (Matrix<int, 12, 2>() << 0b000, 0b001,
                0b000, 0b010,
                0b000, 0b100,
                0b001, 0b011,
                0b001, 0b101,
                0b010, 0b011,
                0b010, 0b110,
                0b011, 0b111,
                0b100, 0b101,
                0b100, 0b110,
                0b101, 0b111,
                0b110, 0b111).finished();

        const Vector3i AXIS_X = Vector3i{1, 0, 0};
        const Vector3i AXIS_Y = Vector3i{0, 1, 0};
        const Vector3i AXIS_Z = Vector3i{0, 0, 1};

        [[nodiscard]]
        Vector3f estimate_centroid(const VectorXf &field_corners) const;

        std::vector<Vector3f> create_vertices(
                MatrixXi &indices,
                const MatrixXf &domain,
                const MatrixXf &sdf,
                int resolution,
                const AlignedBox3f &bounds,
                const NdToFlatIndexer &linearize
        ) const;

        [[nodiscard]]
        std::vector<Vector4i> create_faces(
                const MatrixXi &indices,
                const MatrixXf &sdf,
                int resolution,
                const NdToFlatIndexer &linearize
        ) const;
    };
}