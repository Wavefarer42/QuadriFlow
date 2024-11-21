#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>

#include "entities.h"

using namespace Eigen;

namespace mathext {
    Eigen::MatrixXf clip(
            const Eigen::MatrixXf &mat,
            float minVal,
            float maxVal
    );

    float frobenius_norm_off_diagonal(
            const Eigen::MatrixXf &A
    );

    float percentile(
            const Eigen::VectorXf &vec,
            float percentile
    );

    // Mesh math
    Vector3f face_centroid(
            entities::Mesh &mesh,
            entities::Mesh::VertexFaceIter &face
    );

    MatrixXf face_centroids_ring(
            entities::Mesh &mesh,
            const entities::Mesh::VertexHandle vertex
    );

    MatrixXf count_unique(
            const MatrixXf &mat
    );
}