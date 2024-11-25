#pragma once

#include <Eigen/Core>

#include "entities.h"

using namespace Eigen;

namespace mathext {
    MatrixXf clip(
        const MatrixXf &mat,
        float minVal,
        float maxVal
    );

    float frobenius_norm_off_diagonal(
        const MatrixXf &A
    );

    float percentile(
        const VectorXf &vec,
        float percentile
    );

    // Mesh math
    Vector3f face_centroid(
        entities::Mesh &mesh,
        const entities::Mesh::FaceHandle &face
    );

    MatrixXf face_centroids_ring(
        entities::Mesh &mesh,
        const entities::Mesh::VertexHandle vertex
    );

    MatrixXf count_unique(
        const MatrixXf &mat
    );

    Vector4f intersect_planes(
        const MatrixXf &vertices,
        const MatrixXf &normals
    );
}
