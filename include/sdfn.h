#pragma once

#include <functional>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <tbb/parallel_for.h>

#include "entities.h"

namespace sdfn {
    using namespace Eigen;

    VectorXf sphere(const MatrixXf &domain);

    VectorXf box(const MatrixXf &domain);

    VectorXf cylinder(const MatrixXf &domain);

    entities::SDFn rotate(
            entities::SDFn sdfn,
            const Vector3f axis,
            float angle
    );

    MatrixXf gradient_of(
            const entities::SDFn &sdfn,
            const MatrixXf &domain,
            const float epsilon = 1e-5
    );

    MatrixXf normal_of(
            const entities::SDFn &sdfn,
            const MatrixXf &domain,
            const float epsilon
    );

    MatrixXf normal_of(
            const MatrixXf &gradients
    );
}


