#pragma once

#include <Eigen/Dense>

#include <OpenMesh/Tools/Utils/MeshCheckerT.hh>

namespace svd {
    using namespace Eigen;

    Vector4f intersect_planes(
            const MatrixXf &vertices,
            const MatrixXf &normals
    );
}