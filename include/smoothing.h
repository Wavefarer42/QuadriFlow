#pragma once

#include <Eigen/Core>
#include "entities.h"

namespace smoothing {
    using namespace Eigen;

    // Fields

    VectorXf laplacian_angular_field(
        const entities::SDFn &sdfn,
        entities::Mesh &mesh
    );

    // Smoothing

    entities::Mesh sdfn_projection(
        const entities::SDFn &sdfn,
        entities::Mesh &mesh,
        int iterations = 3,
        float rate = 0.1
    );

    entities::Mesh edge_snapping(
        const entities::SDFn &sdfn,
        entities::Mesh &mesh,
        int iterations = 3,
        float threshold_angle = 30,
        float max_error = 1e-1
    );

    entities::Mesh laplacian_with_sdfn_projection(
        const entities::SDFn &sdfn,
        entities::Mesh &mesh,
        int iterations = 3
    );
};
