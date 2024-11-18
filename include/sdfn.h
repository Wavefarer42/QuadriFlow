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
}


