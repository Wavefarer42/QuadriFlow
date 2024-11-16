#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <tbb/parallel_for.h>

namespace sdfn {
    using namespace Eigen;

    VectorXf sphere(const MatrixXf &domain);

    VectorXf box(const MatrixXf domain);
}


