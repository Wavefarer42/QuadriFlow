#pragma once

#include <vector>

#include <Eigen/Core>

namespace qflow {

    using namespace Eigen;

    void load(const char *filename, MatrixXd &V, MatrixXi &F);

}
