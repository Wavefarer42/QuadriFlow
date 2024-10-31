#pragma once

#include <map>
#include <vector>

#include "entities.h"

#define RCPOVERFLOW 2.93873587705571876e-39f

namespace services {
    using namespace Eigen;

    class Hierarchy {
    public:
        Hierarchy();

        void Initialize(double scale, int with_scale = 0);

        void DownsampleGraph(
                const entities::AdjacentMatrix adj,
                const MatrixXd &V,
                const MatrixXd &N,
                const VectorXd &A,
                MatrixXd &V_p,
                MatrixXd &N_p,
                VectorXd &A_p,
                MatrixXi &to_upper,
                VectorXi &to_lower,
                entities::AdjacentMatrix &adj_p
        );

        void generate_graph_coloring_deterministic(
                const entities::AdjacentMatrix &adj,
                int size,
                std::vector<std::vector<int>> &phases
        );

        void FixFlip();

        void PropagateEdge();

        void DownsampleEdgeGraph(
                std::vector<Vector3i> &FQ,
                std::vector<Vector3i> &F2E,
                std::vector<Vector2i> &edge_diff,
                std::vector<int> &allow_changes,
                int level
        );

        void UpdateGraphValue(
                std::vector<Vector3i> &FQ,
                std::vector<Vector3i> &F2E,
                std::vector<Vector2i> &edge_diff
        );

        enum {
            MAX_DEPTH = 25
        };

        void clearConstraints();

        void propagateConstraints();

        double m_scale;
        int rng_seed;

        MatrixXi m_faces;    // m_faces(i, j) i \in [0, 3) ith index in face j
        VectorXi m_E2E;  // inverse edge
        std::vector<entities::AdjacentMatrix> m_adjacency;
        std::vector<MatrixXd> m_vertices;
        std::vector<MatrixXd> m_normals;
        std::vector<VectorXd> m_vertex_area;
        std::vector<std::vector<std::vector<int>>> m_phases;

        // parameters
        std::vector<MatrixXd> m_orientation;
        std::vector<MatrixXd> m_positions;
        std::vector<VectorXi> mToLower;
        std::vector<MatrixXi> mToUpper;  // mToUpper[h](i, j) \in m_vertices; i \in [0, 2); j \in m_vertices
        std::vector<MatrixXd> m_scales;
        std::vector<MatrixXd> m_areas;

        // constraints
        std::vector<MatrixXd> m_orientation_constraint;
        std::vector<MatrixXd> m_position_constraints;
        std::vector<VectorXd> m_orientation_constraint_weight;
        std::vector<VectorXd> m_position_constraint_weights;

        int with_scale;

        // upper: fine to coarse
        std::vector<std::vector<int>> mToUpperFaces;  // face correspondance
        std::vector<std::vector<int>> mSing;
        std::vector<std::vector<int>> mToUpperEdges; // edge correspondance
        std::vector<std::vector<int>> mToUpperOrients; // rotation of edges from fine to coarse
        std::vector<std::vector<Vector3i>> mFQ; // m_face_edge_orientation
        std::vector<std::vector<Vector3i>> mF2E; // m_face_edge_ids
        std::vector<std::vector<Vector2i>> mE2F; // undirect edges to face ID
        std::vector<std::vector<int> > mAllowChanges;
        std::vector<std::vector<Vector2i>> mEdgeDiff; // face_edgeDiff

    };

}
