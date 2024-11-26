#include <atomic>
#include <fstream>
#include <vector>

#include "dedge.h"

namespace services {
    using namespace Eigen;

    const int INVALID = -1;

    inline int dedge_prev(int e, int deg) { return (e % deg == 0u) ? e + (deg - 1) : e - 1; }

    bool compute_direct_graph(
            MatrixXd &V,
            MatrixXi &F,
            VectorXi &V2E,
            VectorXi &E2E,
            VectorXi &boundary,
            VectorXi &nonManifold
    ) {
        V2E.resize(V.cols());
        V2E.setConstant(INVALID);

        uint32_t deg = F.rows();
        std::vector<std::pair<uint32_t, uint32_t>> tmp(F.size());

        for (int f = 0; f < F.cols(); ++f) {
            for (unsigned int i = 0; i < deg; ++i) {
                unsigned int idx_cur = F(i, f), idx_next = F((i + 1) % deg, f), edge_id = deg * f + i;
                if (idx_cur >= V.cols() || idx_next >= V.cols())
                    throw std::runtime_error("Mesh data contains an out-of-bounds vertex reference!");
                if (idx_cur == idx_next) continue;

                tmp[edge_id] = std::make_pair(idx_next, -1);
                if (V2E[idx_cur] == -1)
                    V2E[idx_cur] = edge_id;
                else {
                    unsigned int idx = V2E[idx_cur];
                    while (tmp[idx].second != -1) {
                        idx = tmp[idx].second;
                    }
                    tmp[idx].second = edge_id;
                }
            }
        }

        nonManifold.resize(V.cols());
        nonManifold.setConstant(false);

        E2E.resize(F.cols() * deg);
        E2E.setConstant(INVALID);

        for (int f = 0; f < F.cols(); ++f) {
            for (uint32_t i = 0; i < deg; ++i) {
                uint32_t idx_cur = F(i, f), idx_next = F((i + 1) % deg, f), edge_id_cur = deg * f + i;

                if (idx_cur == idx_next) continue;

                uint32_t it = V2E[idx_next], edge_id_opp = INVALID;
                while (it != INVALID) {
                    if (tmp[it].first == idx_cur) {
                        if (edge_id_opp == INVALID) {
                            edge_id_opp = it;
                        } else {
                            nonManifold[idx_cur] = true;
                            nonManifold[idx_next] = true;
                            edge_id_opp = INVALID;
                            break;
                        }
                    }
                    it = tmp[it].second;
                }

                if (edge_id_opp != INVALID && edge_id_cur < edge_id_opp) {
                    E2E[edge_id_cur] = edge_id_opp;
                    E2E[edge_id_opp] = edge_id_cur;
                }
            }
        }
        std::atomic<uint32_t> nonManifoldCounter(0), boundaryCounter(0), isolatedCounter(0);

        boundary.resize(V.cols());
        boundary.setConstant(false);

        /* Detect boundary regions of the mesh and adjust vertex->edge pointers*/
        for (int i = 0; i < V.cols(); ++i) {
            uint32_t edge = V2E[i];
            if (edge == INVALID) {
                isolatedCounter++;
                continue;
            }
            if (nonManifold[i]) {
                nonManifoldCounter++;
                V2E[i] = INVALID;
                continue;
            }

            /* Walk backwards to the first boundary edge (if any) */
            uint32_t start = edge, v2e = INVALID;
            do {
                v2e = std::min(v2e, edge);
                uint32_t prevEdge = E2E[dedge_prev(edge, deg)];
                if (prevEdge == INVALID) {
                    /* Reached boundary -- update the vertex->edge link */
                    v2e = edge;
                    boundary[i] = true;
                    boundaryCounter++;
                    break;
                }
                edge = prevEdge;
            } while (edge != start);
            V2E[i] = v2e;
        }

        return true;
    }

    void compute_direct_graph_quad(
            std::vector<Vector3d> &V,
            std::vector<Vector4i> &F,
            std::vector<int> &V2E,
            std::vector<int> &E2E,
            VectorXi &boundary,
            VectorXi &nonManifold
    ) {
        V2E.clear();
        E2E.clear();
        boundary = VectorXi();
        nonManifold = VectorXi();
        V2E.resize(V.size(), INVALID);

        uint32_t deg = 4;
        std::vector<std::pair<uint32_t, uint32_t>> tmp(F.size() * deg);

        for (int f = 0; f < F.size(); ++f) {
            for (unsigned int i = 0; i < deg; ++i) {
                unsigned int idx_cur = F[f][i], idx_next = F[f][(i + 1) % deg], edge_id = deg * f + i;
                if (idx_cur >= V.size() || idx_next >= V.size())
                    throw std::runtime_error("Mesh data contains an out-of-bounds vertex reference!");
                if (idx_cur == idx_next) continue;
                tmp[edge_id] = std::make_pair(idx_next, -1);
                if (V2E[idx_cur] == -1) {
                    V2E[idx_cur] = edge_id;
                } else {
                    unsigned int idx = V2E[idx_cur];
                    while (tmp[idx].second != -1) {
                        idx = tmp[idx].second;
                    }
                    tmp[idx].second = edge_id;
                }
            }
        }

        nonManifold.resize(V.size());
        nonManifold.setConstant(false);

        E2E.resize(F.size() * deg, INVALID);

        for (int f = 0; f < F.size(); ++f) {
            for (uint32_t i = 0; i < deg; ++i) {
                uint32_t idx_cur = F[f][i], idx_next = F[f][(i + 1) % deg], edge_id_cur = deg * f + i;

                if (idx_cur == idx_next) continue;

                uint32_t it = V2E[idx_next], edge_id_opp = INVALID;
                while (it != INVALID) {
                    if (tmp[it].first == idx_cur) {
                        if (edge_id_opp == INVALID) {
                            edge_id_opp = it;
                        } else {
                            nonManifold[idx_cur] = true;
                            nonManifold[idx_next] = true;
                            edge_id_opp = INVALID;
                            break;
                        }
                    }
                    it = tmp[it].second;
                }

                if (edge_id_opp != INVALID && edge_id_cur < edge_id_opp) {
                    E2E[edge_id_cur] = edge_id_opp;
                    E2E[edge_id_opp] = edge_id_cur;
                }
            }
        }
        std::atomic<uint32_t> nonManifoldCounter(0), boundaryCounter(0), isolatedCounter(0);

        boundary.resize(V.size());
        boundary.setConstant(false);

        /* Detect boundary regions of the mesh and adjust vertex->edge pointers*/
        for (int i = 0; i < V.size(); ++i) {
            uint32_t edge = V2E[i];
            if (edge == INVALID) {
                isolatedCounter++;
                continue;
            }
            if (nonManifold[i]) {
                nonManifoldCounter++;
                V2E[i] = INVALID;
                continue;
            }

            /* Walk backwards to the first boundary edge (if any) */
            uint32_t start = edge, v2e = INVALID;
            do {
                v2e = std::min(v2e, edge);
                uint32_t prevEdge = E2E[dedge_prev(edge, deg)];
                if (prevEdge == INVALID) {
                    /* Reached boundary -- update the vertex->edge link */
                    v2e = edge;
                    boundary[i] = true;
                    boundaryCounter++;
                    break;
                }
                edge = prevEdge;
            } while (edge != start);
            V2E[i] = v2e;
        }
//        printf("counter %d %d\n", (int)boundaryCounter, (int)nonManifoldCounter);
    }
}
