#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include "config.hpp"
#include "field-math.hpp"
#include "hierarchy.hpp"

namespace qflow {

    class Optimizer {
    public:
        Optimizer();

        static void optimize_orientations(Hierarchy &mRes);

        static void optimize_scale(Hierarchy &mRes, VectorXd &rho, int adaptive);

        static void optimize_positions(Hierarchy &mRes, int with_scale = 0);

        static void optimize_integer_constraints(Hierarchy &mRes, std::map<int, int> &singularities);

        static void optimize_positions_fixed(
                Hierarchy &mRes, std::vector<DEdge> &edge_values, std::vector<Vector2i> &edge_diff,
                std::set<int> &sharp_vertices,
                std::map<int, std::pair<Vector3d, Vector3d>> &sharp_constraints, int with_scale = 0);

        static void optimize_positions_sharp(
                Hierarchy &mRes, std::vector<DEdge> &edge_values, std::vector<Vector2i> &edge_diff,
                std::vector<int> &sharp_edges, std::set<int> &sharp_vertices,
                std::map<int, std::pair<Vector3d, Vector3d>> &sharp_constraints, int with_scale = 0);

        static void optimize_positions_dynamic(
                MatrixXi &F, MatrixXd &V, MatrixXd &N, MatrixXd &Q, std::vector<std::vector<int>> &Vset,
                std::vector<Vector3d> &O_compact, std::vector<Vector4i> &F_compact,
                std::vector<int> &V2E_compact, std::vector<int> &E2E_compact, double mScale,
                std::vector<Vector3d> &diffs, std::vector<int> &diff_count,
                std::map<std::pair<int, int>, int> &o2e, std::vector<int> &sharp_o,
                std::map<int, std::pair<Vector3d, Vector3d>> &compact_sharp_constraints, int with_scale);
    };

} // namespace qflow

#endif
