#pragma once

#include <Eigen/Core>

using namespace Eigen;

namespace transformations {
    Vector4f intersect_planes(
            const MatrixXf &vertices,
            const MatrixXf &normals
    );
}