#pragma once

#include "entities.h"

using namespace Eigen;

namespace surfacenets {

    MatrixXf cartesian_product(
            const VectorXf &X,
            const VectorXf &Y,
            const VectorXf &Z
    );

    MatrixXf create_sample_grid(
            const AlignedBox3f &bounds,
            float resolution
    );

    AlignedBox3f estimate_bounding_box(
            entities::SDFn sdfn,
            int grid_resolution
    );

    class SurfaceNetsMeshStrategy {
    public:
        entities::QuadMesh mesh(entities::SDFn sdfn, int resolution) const;

    private:
        const Matrix<int, 8, 3> CUBE_CORNERS = (Matrix<int, 8, 3>() << 0, 0, 0,
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

        const Matrix<int, 3, 3> AXIS = (Matrix<int, 3, 3>() << 1, 0, 0,
                0, 1, 0,
                0, 0, 1).finished();
    };
}