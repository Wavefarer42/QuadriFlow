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
        const MatrixXf &normals,
        int patience = 10
    );

    /**
    * Takes a column-major matrix of vertices and normalizes them between -1 and 1.
    */
    std::tuple<MatrixXf, float, Vector3f> normalize(
        const MatrixXf &vertices
    );

    /**
     * Denormalizes the vertices with the given scale and offset.
     */
    MatrixXf denormalize(
        const MatrixXf &vertices,
        const float scale,
        const Vector3f &offset
    );
}
