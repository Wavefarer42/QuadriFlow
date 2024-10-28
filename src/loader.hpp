#ifndef __LOADER_H
#define __LOADER_H

#include <vector>

#include <Eigen/Core>

namespace qflow {

    using namespace Eigen;

    void load(const char *filename, MatrixXd &V, MatrixXi &F);

} // namespace qflow

#endif
