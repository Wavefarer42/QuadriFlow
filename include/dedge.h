#pragma once

#include <vector>

#include "entities.h"

namespace services {
    using namespace Eigen;

    inline int dedge_prev_3(int e) { return (e % 3 == 0) ? e + 2 : e - 1; }

    inline int dedge_next_3(int e) { return (e % 3 == 2) ? e - 2 : e + 1; }

    bool compute_direct_graph(
            MatrixXd &V,
            MatrixXi &F,
            VectorXi &V2E,
            VectorXi &E2E,
            VectorXi &boundary,
            VectorXi &nonManifold
    );

    void compute_direct_graph_quad(
            std::vector<Vector3d> &V,
            std::vector<Vector4i> &F,
            std::vector<int> &V2E,
            std::vector<int> &E2E,
            VectorXi &boundary,
            VectorXi &nonManifold
    );

}
