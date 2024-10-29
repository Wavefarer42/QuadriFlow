#pragma once

#include <vector>

#include <Eigen/Core>

namespace qflow {

    using namespace Eigen;

    enum class SolverStatus {
        Sat,
        Unsat,
        Timeout,
    };

    SolverStatus SolveSatProblem(
            int n_variable,
            std::vector<int> &value,
            const std::vector<bool> flexible,  // NOQA
            const std::vector<Vector3i> &variable_eq,
            const std::vector<Vector3i> &constant_eq,
            const std::vector<Vector4i> &variable_ge,
            const std::vector<Vector2i> &constant_ge,
            int timeout = 8
    );
}
