#include <algorithm>
#include <queue>
#include <unordered_map>
#include <cmath>
#include <fstream>
#include <memory>
#include <Eigen/Sparse>
#include <list>
#include <map>
#include <set>
#include <random>
#include <mach/mig.h>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include "pcg32.h"

#include "quadriflow.h"
#include "entities.h"
#include "smoothing.h"
#include "mathext.h"
#include "sdfn.h"


using namespace Eigen;


namespace quadriflow {
    Hierarchy::Hierarchy() {
        m_adjacency.resize(MAX_DEPTH + 1);
        m_vertices.resize(MAX_DEPTH + 1);
        m_normals.resize(MAX_DEPTH + 1);
        m_vertex_area.resize(MAX_DEPTH + 1);
        m_phases.resize(MAX_DEPTH + 1);
        mToLower.resize(MAX_DEPTH);
        mToUpper.resize(MAX_DEPTH);
        rng_seed = 0;

        m_orientation_constraint.reserve(MAX_DEPTH + 1);
        m_orientation_constraint_weight.reserve(MAX_DEPTH + 1);
        m_position_constraints.reserve(MAX_DEPTH + 1);
        m_position_constraint_weights.reserve(MAX_DEPTH + 1);
    }

    void Hierarchy::Initialize(double scale, int with_scale) {
        this->with_scale = with_scale;
        generate_graph_coloring_deterministic(m_adjacency[0], m_vertices[0].cols(), m_phases[0]);

        for (int i = 0; i < MAX_DEPTH; ++i) {
            DownsampleGraph(m_adjacency[i], m_vertices[i], m_normals[i], m_vertex_area[i], m_vertices[i + 1],
                            m_normals[i + 1], m_vertex_area[i + 1], mToUpper[i],
                            mToLower[i], m_adjacency[i + 1]);
            generate_graph_coloring_deterministic(m_adjacency[i + 1], m_vertices[i + 1].cols(), m_phases[i + 1]);
            if (m_vertices[i + 1].cols() == 1) {
                m_adjacency.resize(i + 2);
                m_vertices.resize(i + 2);
                m_normals.resize(i + 2);
                m_vertex_area.resize(i + 2);
                mToUpper.resize(i + 1);
                mToLower.resize(i + 1);
                break;
            }
        }
        m_orientation.resize(m_vertices.size());
        m_positions.resize(m_vertices.size());
        m_scales.resize(m_vertices.size());
        m_areas.resize(m_vertices.size());

        m_position_constraints.resize(m_vertices.size());
        m_position_constraint_weights.resize(m_vertices.size());
        m_orientation_constraint.resize(m_vertices.size());
        m_orientation_constraint_weight.resize(m_vertices.size());

        // Set random seed
        srand(rng_seed);

        m_scale = scale;
        for (int i = 0; i < m_vertices.size(); ++i) {
            m_orientation[i].resize(m_normals[i].rows(), m_normals[i].cols());
            m_positions[i].resize(m_normals[i].rows(), m_normals[i].cols());
            m_scales[i].resize(2, m_normals[i].cols());
            m_areas[i].resize(2, m_normals[i].cols());
            for (int j = 0; j < m_normals[i].cols(); ++j) {
                Vector3d s, t;
                mathext::coordinate_system(m_normals[i].col(j), s, t);
                // rand() is not thread safe!
                double angle = ((double) rand()) / RAND_MAX * 2 * M_PI;
                double x = ((double) rand()) / RAND_MAX * 2 - 1.f;
                double y = ((double) rand()) / RAND_MAX * 2 - 1.f;
                m_orientation[i].col(j) = s * std::cos(angle) + t * std::sin(angle);
                m_positions[i].col(j) = m_vertices[i].col(j) + (s * x + t * y) * scale;
                if (with_scale) {
                    m_scales[i].col(j) = Vector2d(1.0f, 1.0f);
                    m_areas[i].col(j) = Vector2d(0.0, 0.0);
                }
            }
        }
    }

    void Hierarchy::generate_graph_coloring_deterministic(
        const AdjacentMatrix &adj,
        int size,
        std::vector<std::vector<int> > &phases
    ) {
        phases.clear();

        std::vector<uint32_t> perm(size);
        for (uint32_t i = 0; i < size; ++i) perm[i] = i;
        pcg32 rng;
        rng.shuffle(perm.begin(), perm.end());

        std::vector<int> color(size, -1);
        std::vector<uint8_t> possible_colors;
        std::vector<int> size_per_color;
        int ncolors = 0;

        for (uint32_t i = 0; i < size; ++i) {
            uint32_t ip = perm[i];

            std::fill(possible_colors.begin(), possible_colors.end(), 1);

            for (auto &link: adj[ip]) {
                int c = color[link.id];
                if (c >= 0) possible_colors[c] = 0;
            }

            int chosen_color = -1;
            for (uint32_t j = 0; j < possible_colors.size(); ++j) {
                if (possible_colors[j]) {
                    chosen_color = j;
                    break;
                }
            }

            if (chosen_color < 0) {
                chosen_color = ncolors++;
                possible_colors.resize(ncolors);
                size_per_color.push_back(0);
            }

            color[ip] = chosen_color;
            size_per_color[chosen_color]++;
        }
        phases.resize(ncolors);
        for (int i = 0; i < ncolors; ++i) phases[i].reserve(size_per_color[i]);
        for (uint32_t i = 0; i < size; ++i) phases[color[i]].push_back(i);
    }

    void Hierarchy::DownsampleGraph(
        const AdjacentMatrix adj,
        const MatrixXd &V,
        const MatrixXd &N,
        const VectorXd &A,
        MatrixXd &V_p,
        MatrixXd &N_p,
        VectorXd &A_p,
        MatrixXi &to_upper,
        VectorXi &to_lower,
        AdjacentMatrix &adj_p
    ) {
        struct Entry {
            int i, j;
            double order;

            inline Entry() { i = j = -1; };

            inline Entry(int i, int j, double order) : i(i), j(j), order(order) {
            }

            inline bool operator<(const Entry &e) const { return order > e.order; }

            inline bool operator==(const Entry &e) const { return order == e.order; }
        };

        int nLinks = 0;
        for (auto &adj_i: adj) nLinks += adj_i.size();
        std::vector<Entry> entries(nLinks);
        std::vector<int> bases(adj.size());
        for (int i = 1; i < bases.size(); ++i) {
            bases[i] = bases[i - 1] + adj[i - 1].size();
        }

        for (int i = 0; i < V.cols(); ++i) {
            int base = bases[i];
            auto &ad = adj[i];
            auto entry_it = entries.begin() + base;
            for (auto it = ad.begin(); it != ad.end(); ++it, ++entry_it) {
                int k = it->id;
                double dp = N.col(i).dot(N.col(k));
                double ratio = A[i] > A[k] ? (A[i] / A[k]) : (A[k] / A[i]);
                *entry_it = Entry(i, k, dp * ratio);
            }
        }

        std::stable_sort(entries.begin(), entries.end(), std::less<Entry>());

        std::vector<bool> mergeFlag(V.cols(), false);

        int nCollapsed = 0;
        for (int i = 0; i < nLinks; ++i) {
            const Entry &e = entries[i];
            if (mergeFlag[e.i] || mergeFlag[e.j]) continue;
            mergeFlag[e.i] = mergeFlag[e.j] = true;
            entries[nCollapsed++] = entries[i];
        }
        int vertexCount = V.cols() - nCollapsed;

        // Allocate memory for coarsened graph
        V_p.resize(3, vertexCount);
        N_p.resize(3, vertexCount);
        A_p.resize(vertexCount);
        to_upper.resize(2, vertexCount);
        to_lower.resize(V.cols());

        for (int i = 0; i < nCollapsed; ++i) {
            const Entry &e = entries[i];
            const double area1 = A[e.i], area2 = A[e.j], surfaceArea = area1 + area2;
            if (surfaceArea > RCPOVERFLOW)
                V_p.col(i) = (V.col(e.i) * area1 + V.col(e.j) * area2) / surfaceArea;
            else
                V_p.col(i) = (V.col(e.i) + V.col(e.j)) * 0.5f;
            Vector3d normal = N.col(e.i) * area1 + N.col(e.j) * area2;
            double norm = normal.norm();
            N_p.col(i) = norm > RCPOVERFLOW ? Vector3d(normal / norm) : Vector3d::UnitX();
            A_p[i] = surfaceArea;
            to_upper.col(i) << e.i, e.j;
            to_lower[e.i] = i;
            to_lower[e.j] = i;
        }

        int offset = nCollapsed;

        for (int i = 0; i < V.cols(); ++i) {
            if (!mergeFlag[i]) {
                int idx = offset++;
                V_p.col(idx) = V.col(i);
                N_p.col(idx) = N.col(i);
                A_p[idx] = A[i];
                to_upper.col(idx) << i, -1;
                to_lower[i] = idx;
            }
        }

        adj_p.resize(V_p.cols());
        std::vector<int> capacity(V_p.cols());
        std::vector<std::vector<Link> > scratches(V_p.cols());
        for (int i = 0; i < V_p.cols(); ++i) {
            int t = 0;
            for (int j = 0; j < 2; ++j) {
                int upper = to_upper(j, i);
                if (upper == -1) continue;
                t += adj[upper].size();
            }
            scratches[i].reserve(t);
            adj_p[i].reserve(t);
        }
        for (int i = 0; i < V_p.cols(); ++i) {
            auto &scratch = scratches[i];
            for (int j = 0; j < 2; ++j) {
                int upper = to_upper(j, i);
                if (upper == -1) continue;
                auto &ad = adj[upper];
                for (auto &link: ad) scratch.push_back(Link(to_lower[link.id], link.weight));
            }
            std::sort(scratch.begin(), scratch.end());
            int id = -1;
            auto &ad = adj_p[i];
            for (auto &link: scratch) {
                if (link.id != i) {
                    if (id != link.id) {
                        ad.push_back(link);
                        id = link.id;
                    } else {
                        ad.back().weight += link.weight;
                    }
                }
            }
        }
    }

    void Hierarchy::UpdateGraphValue(
        std::vector<Vector3i> &FQ,
        std::vector<Vector3i> &F2E,
        std::vector<Vector2i> &edge_diff
    ) {
        FQ = std::move(mFQ[0]);
        F2E = std::move(mF2E[0]);
        edge_diff = std::move(mEdgeDiff[0]);
    }

    void Hierarchy::DownsampleEdgeGraph(
        std::vector<Vector3i> &FQ,
        std::vector<Vector3i> &F2E,
        std::vector<Vector2i> &edge_diff,
        std::vector<int> &allow_changes,
        int level
    ) {
        std::vector<Vector2i> E2F(edge_diff.size(), Vector2i(-1, -1));
        for (int i = 0; i < F2E.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                int e = F2E[i][j];
                if (E2F[e][0] == -1)
                    E2F[e][0] = i;
                else
                    E2F[e][1] = i;
            }
        }
        int levels = (level == -1) ? 100 : level;
        mFQ.resize(levels);
        mF2E.resize(levels);
        mE2F.resize(levels);
        mEdgeDiff.resize(levels);
        mAllowChanges.resize(levels);
        mSing.resize(levels);
        mToUpperEdges.resize(levels - 1);
        mToUpperOrients.resize(levels - 1);
        for (int i = 0; i < FQ.size(); ++i) {
            Vector2i diff(0, 0);
            for (int j = 0; j < 3; ++j) {
                diff += mathext::rshift90(edge_diff[F2E[i][j]], FQ[i][j]);
            }
            if (diff != Vector2i::Zero()) {
                mSing[0].push_back(i);
            }
        }
        mAllowChanges[0] = allow_changes;
        mFQ[0] = std::move(FQ);
        mF2E[0] = std::move(F2E);
        mE2F[0] = std::move(E2F);
        mEdgeDiff[0] = std::move(edge_diff);
        for (int l = 0; l < levels - 1; ++l) {
            auto &FQ = mFQ[l];
            auto &E2F = mE2F[l];
            auto &F2E = mF2E[l];
            auto &Allow = mAllowChanges[l];
            auto &EdgeDiff = mEdgeDiff[l];
            auto &Sing = mSing[l];
            std::vector<int> fixed_faces(F2E.size(), 0);
            for (auto &s: Sing) {
                fixed_faces[s] = 1;
            }

            auto &toUpper = mToUpperEdges[l];
            auto &toUpperOrients = mToUpperOrients[l];
            toUpper.resize(E2F.size(), -1);
            toUpperOrients.resize(E2F.size(), 0);

            auto &nFQ = mFQ[l + 1];
            auto &nE2F = mE2F[l + 1];
            auto &nF2E = mF2E[l + 1];
            auto &nAllow = mAllowChanges[l + 1];
            auto &nEdgeDiff = mEdgeDiff[l + 1];
            auto &nSing = mSing[l + 1];

            for (int i = 0; i < E2F.size(); ++i) {
                if (EdgeDiff[i] != Vector2i::Zero()) continue;
                if ((E2F[i][0] >= 0 && fixed_faces[E2F[i][0]]) ||
                    (E2F[i][1] >= 0 && fixed_faces[E2F[i][1]])) {
                    continue;
                }
                for (int j = 0; j < 2; ++j) {
                    int f = E2F[i][j];
                    if (f < 0) continue;
                    for (int k = 0; k < 3; ++k) {
                        int neighbor_e = F2E[f][k];
                        for (int m = 0; m < 2; ++m) {
                            int neighbor_f = E2F[neighbor_e][m];
                            if (neighbor_f < 0) continue;
                            if (fixed_faces[neighbor_f] == 0) fixed_faces[neighbor_f] = 1;
                        }
                    }
                }
                if (E2F[i][0] >= 0) fixed_faces[E2F[i][0]] = 2;
                if (E2F[i][1] >= 0) fixed_faces[E2F[i][1]] = 2;
                toUpper[i] = -2;
            }
            for (int i = 0; i < E2F.size(); ++i) {
                if (toUpper[i] == -2) continue;
                if ((E2F[i][0] < 0 || fixed_faces[E2F[i][0]] == 2) &&
                    (E2F[i][1] < 0 || fixed_faces[E2F[i][1]] == 2)) {
                    toUpper[i] = -3;
                    continue;
                }
            }
            int numE = 0;
            for (int i = 0; i < toUpper.size(); ++i) {
                if (toUpper[i] == -1) {
                    if ((E2F[i][0] < 0 || fixed_faces[E2F[i][0]] < 2) &&
                        (E2F[i][1] < 0 || fixed_faces[E2F[i][1]] < 2)) {
                        nE2F.push_back(E2F[i]);
                        toUpperOrients[i] = 0;
                        toUpper[i] = numE++;
                        continue;
                    }
                    int f0 = (E2F[i][1] < 0 || fixed_faces[E2F[i][0]] < 2) ? E2F[i][0] : E2F[i][1];
                    int e = i;
                    int f = f0;
                    std::vector<std::pair<int, int> > paths;
                    paths.push_back(std::make_pair(i, 0));
                    while (true) {
                        if (E2F[e][0] == f)
                            f = E2F[e][1];
                        else if (E2F[e][1] == f)
                            f = E2F[e][0];
                        if (f < 0 || fixed_faces[f] < 2) {
                            for (int j = 0; j < paths.size(); ++j) {
                                auto &p = paths[j];
                                toUpper[p.first] = numE;
                                int orient = p.second;
                                if (j > 0) orient = (orient + toUpperOrients[paths[j - 1].first]) % 4;
                                toUpperOrients[p.first] = orient;
                            }
                            nE2F.push_back(Vector2i(f0, f));
                            numE += 1;
                            break;
                        }
                        int ind0 = -1, ind1 = -1;
                        int e0 = e;
                        for (int j = 0; j < 3; ++j) {
                            if (F2E[f][j] == e) {
                                ind0 = j;
                                break;
                            }
                        }
                        for (int j = 0; j < 3; ++j) {
                            int e1 = F2E[f][j];
                            if (e1 != e && toUpper[e1] != -2) {
                                e = e1;
                                ind1 = j;
                                break;
                            }
                        }

                        if (ind1 != -1) {
                            paths.push_back(std::make_pair(e, (FQ[f][ind1] - FQ[f][ind0] + 6) % 4));
                        } else {
                            if (EdgeDiff[e] != Vector2i::Zero()) {
                                printf("Unsatisfied !!!...\n");
                                printf("%d %d %d: %d %d\n", F2E[f][0], F2E[f][1], F2E[f][2], e0, e);
                                exit(0);
                            }
                            for (auto &p: paths) {
                                toUpper[p.first] = numE;
                                toUpperOrients[p.first] = 0;
                            }
                            numE += 1;
                            nE2F.push_back(Vector2i(f0, f0));
                            break;
                        }
                    }
                }
            }
            nEdgeDiff.resize(numE);
            nAllow.resize(numE * 2, 1);
            for (int i = 0; i < toUpper.size(); ++i) {
                if (toUpper[i] >= 0 && toUpperOrients[i] == 0) {
                    nEdgeDiff[toUpper[i]] = EdgeDiff[i];
                }
                if (toUpper[i] >= 0) {
                    int dimension = toUpperOrients[i] % 2;
                    if (Allow[i * 2 + dimension] == 0)
                        nAllow[toUpper[i] * 2] = 0;
                    else if (Allow[i * 2 + dimension] == 2)
                        nAllow[toUpper[i] * 2] = 2;
                    if (Allow[i * 2 + 1 - dimension] == 0)
                        nAllow[toUpper[i] * 2 + 1] = 0;
                    else if (Allow[i * 2 + 1 - dimension] == 2)
                        nAllow[toUpper[i] * 2 + 1] = 2;
                }
            }
            std::vector<int> upperface(F2E.size(), -1);

            for (int i = 0; i < F2E.size(); ++i) {
                Vector3i eid;
                for (int j = 0; j < 3; ++j) {
                    eid[j] = toUpper[F2E[i][j]];
                }
                if (eid[0] >= 0 && eid[1] >= 0 && eid[2] >= 0) {
                    Vector3i eid_orient;
                    for (int j = 0; j < 3; ++j) {
                        eid_orient[j] = (FQ[i][j] + 4 - toUpperOrients[F2E[i][j]]) % 4;
                    }
                    upperface[i] = nF2E.size();
                    nF2E.push_back(eid);
                    nFQ.push_back(eid_orient);
                }
            }
            for (int i = 0; i < nE2F.size(); ++i) {
                for (int j = 0; j < 2; ++j) {
                    if (nE2F[i][j] >= 0) nE2F[i][j] = upperface[nE2F[i][j]];
                }
            }

            for (auto &s: Sing) {
                if (upperface[s] >= 0) nSing.push_back(upperface[s]);
            }
            mToUpperFaces.push_back(std::move(upperface));

            if (nEdgeDiff.size() == EdgeDiff.size()) {
                levels = l + 1;
                break;
            }
        }

        mFQ.resize(levels);
        mF2E.resize(levels);
        mAllowChanges.resize(levels);
        mE2F.resize(levels);
        mEdgeDiff.resize(levels);
        mSing.resize(levels);
        mToUpperEdges.resize(levels - 1);
        mToUpperOrients.resize(levels - 1);
    }

    void Hierarchy::FixFlip() {
        int l = mF2E.size() - 1;
        auto &F2E = mF2E[l];
        auto &E2F = mE2F[l];
        auto &FQ = mFQ[l];
        auto &EdgeDiff = mEdgeDiff[l];
        auto &AllowChange = mAllowChanges[l];

        // build m_E2E
        std::vector<int> E2E(F2E.size() * 3, -1);
        for (int i = 0; i < E2F.size(); ++i) {
            int v1 = E2F[i][0];
            int v2 = E2F[i][1];
            int t1 = 0;
            int t2 = 2;
            if (v1 != -1)
                while (F2E[v1][t1] != i) t1 += 1;
            if (v2 != -1)
                while (F2E[v2][t2] != i) t2 -= 1;
            t1 += v1 * 3;
            t2 += v2 * 3;
            if (v1 != -1) E2E[t1] = (v2 == -1) ? -1 : t2;
            if (v2 != -1) E2E[t2] = (v1 == -1) ? -1 : t1;
        }

        auto Area = [&](int f) {
            Vector2i diff1 = mathext::rshift90(EdgeDiff[F2E[f][0]], FQ[f][0]);
            Vector2i diff2 = mathext::rshift90(EdgeDiff[F2E[f][1]], FQ[f][1]);
            return diff1[0] * diff2[1] - diff1[1] * diff2[0];
        };
        std::vector<int> valences(F2E.size() * 3, -10000); // comment this line
        auto CheckShrink = [&](int deid, int allowed_edge_length) {
            // Check if we want shrink direct edge deid so that all edge length is smaller than
            // allowed_edge_length
            if (deid == -1) {
                return false;
            }
            std::vector<int> corresponding_faces;
            std::vector<int> corresponding_edges;
            std::vector<Vector2i> corresponding_diff;
            int deid0 = deid;
            while (deid != -1) {
                deid = deid / 3 * 3 + (deid + 2) % 3;
                if (E2E[deid] == -1) break;
                deid = E2E[deid];
                if (deid == deid0) break;
            }
            Vector2i diff = EdgeDiff[F2E[deid / 3][deid % 3]];
            do {
                corresponding_diff.push_back(diff);
                corresponding_edges.push_back(deid);
                corresponding_faces.push_back(deid / 3);

                // transform to the next face
                deid = E2E[deid];
                if (deid == -1) {
                    return false;
                }
                // transform for the target incremental diff
                diff = -mathext::rshift90(diff, FQ[deid / 3][deid % 3]);
                deid = deid / 3 * 3 + (deid + 1) % 3;
                // transform to local
                diff = mathext::rshift90(diff, (4 - FQ[deid / 3][deid % 3]) % 4);
            } while (deid != corresponding_edges.front());
            // check diff
            if (deid != -1 && diff != corresponding_diff.front()) {
                return false;
            }
            std::unordered_map<int, Vector2i> new_values;
            for (int i = 0; i < corresponding_diff.size(); ++i) {
                int deid = corresponding_edges[i];
                int eid = F2E[deid / 3][deid % 3];
                new_values[eid] = EdgeDiff[eid];
            }
            for (int i = 0; i < corresponding_diff.size(); ++i) {
                int deid = corresponding_edges[i];
                int eid = F2E[deid / 3][deid % 3];
                for (int j = 0; j < 2; ++j) {
                    if (corresponding_diff[i][j] != 0 && AllowChange[eid * 2 + j] == 0) return false;
                }
                auto &res = new_values[eid];
                res -= corresponding_diff[i];
                int edge_thres = allowed_edge_length;
                if (abs(res[0]) > edge_thres || abs(res[1]) > edge_thres) {
                    return false;
                }
                if ((abs(res[0]) > 1 && abs(res[1]) != 0) || (abs(res[1]) > 1 && abs(res[0]) != 0))
                    return false;
            }
            int prev_area = 0, current_area = 0;
            for (int f = 0; f < corresponding_faces.size(); ++f) {
                int area = Area(corresponding_faces[f]);
                if (area < 0) prev_area += 1;
            }
            for (auto &p: new_values) {
                std::swap(EdgeDiff[p.first], p.second);
            }
            for (int f = 0; f < corresponding_faces.size(); ++f) {
                int area = Area(corresponding_faces[f]);
                if (area < 0) {
                    current_area += 1;
                }
            }
            if (current_area < prev_area) {
                return true;
            }
            for (auto &p: new_values) {
                std::swap(EdgeDiff[p.first], p.second);
            }
            return false;
        };

        std::queue<int> flipped;
        for (int i = 0; i < F2E.size(); ++i) {
            int area = Area(i);
            if (area < 0) {
                flipped.push(i);
            }
        }

        bool update = false;
        int max_len = 1;
        while (!update && max_len <= 2) {
            while (!flipped.empty()) {
                int f = flipped.front();
                if (Area(f) >= 0) {
                    flipped.pop();
                    continue;
                }
                for (int i = 0; i < 3; ++i) {
                    if (CheckShrink(f * 3 + i, max_len) || CheckShrink(E2E[f * 3 + i], max_len)) {
                        update = true;
                        break;
                    }
                }
                flipped.pop();
            }
            max_len += 1;
        }
        if (update) {
            Hierarchy flip_hierarchy;
            flip_hierarchy.DownsampleEdgeGraph(mFQ.back(), mF2E.back(), mEdgeDiff.back(),
                                               mAllowChanges.back(), -1);
            flip_hierarchy.FixFlip();
            flip_hierarchy.UpdateGraphValue(mFQ.back(), mF2E.back(), mEdgeDiff.back());
        }
        PropagateEdge();
    }

    void Hierarchy::PropagateEdge() {
        for (int level = mToUpperEdges.size(); level > 0; --level) {
            auto &EdgeDiff = mEdgeDiff[level];
            auto &nEdgeDiff = mEdgeDiff[level - 1];
            auto &FQ = mFQ[level];
            auto &nFQ = mFQ[level - 1];
            auto &F2E = mF2E[level - 1];
            auto &toUpper = mToUpperEdges[level - 1];
            auto &toUpperFace = mToUpperFaces[level - 1];
            auto &toUpperOrients = mToUpperOrients[level - 1];
            for (int i = 0; i < toUpper.size(); ++i) {
                if (toUpper[i] >= 0) {
                    int orient = (4 - toUpperOrients[i]) % 4;
                    nEdgeDiff[i] = mathext::rshift90(EdgeDiff[toUpper[i]], orient);
                } else {
                    nEdgeDiff[i] = Vector2i(0, 0);
                }
            }
            for (int i = 0; i < toUpperFace.size(); ++i) {
                if (toUpperFace[i] == -1) continue;
                Vector3i eid_orient = FQ[toUpperFace[i]];
                for (int j = 0; j < 3; ++j) {
                    nFQ[i][j] = (eid_orient[j] + toUpperOrients[F2E[i][j]]) % 4;
                }
            }
        }
    }

    void Hierarchy::clearConstraints() {
        int levels = m_vertices.size();
        if (levels == 0) return;
        for (int i = 0; i < levels; ++i) {
            int size = m_vertices[i].cols();
            m_orientation_constraint[i].resize(3, size);
            m_position_constraints[i].resize(3, size);
            m_orientation_constraint_weight[i].resize(size);
            m_position_constraint_weights[i].resize(size);
            m_orientation_constraint_weight[i].setZero();
            m_position_constraint_weights[i].setZero();
        }
    }

    void Hierarchy::propagateConstraints() {
        int levels = m_vertices.size();
        if (levels == 0) return;

        for (int l = 0; l < levels - 1; ++l) {
            auto &N = m_normals[l];
            auto &N_next = m_normals[l + 1];
            auto &V = m_vertices[l];
            auto &V_next = m_vertices[l + 1];
            auto &CQ = m_orientation_constraint[l];
            auto &CQ_next = m_orientation_constraint[l + 1];
            auto &CQw = m_orientation_constraint_weight[l];
            auto &CQw_next = m_orientation_constraint_weight[l + 1];
            auto &CO = m_position_constraints[l];
            auto &CO_next = m_position_constraints[l + 1];
            auto &COw = m_position_constraint_weights[l];
            auto &COw_next = m_position_constraint_weights[l + 1];
            auto &toUpper = mToUpper[l];
            MatrixXd &S = m_scales[l];

            for (uint32_t i = 0; i != m_vertices[l + 1].cols(); ++i) {
                Vector2i upper = toUpper.col(i);
                Vector3d cq = Vector3d::Zero(), co = Vector3d::Zero();
                float cqw = 0.0f, cow = 0.0f;

                bool has_cq0 = CQw[upper[0]] != 0;
                bool has_cq1 = upper[1] != -1 && CQw[upper[1]] != 0;
                bool has_co0 = COw[upper[0]] != 0;
                bool has_co1 = upper[1] != -1 && COw[upper[1]] != 0;

                if (has_cq0 && !has_cq1) {
                    cq = CQ.col(upper[0]);
                    cqw = CQw[upper[0]];
                } else if (has_cq1 && !has_cq0) {
                    cq = CQ.col(upper[1]);
                    cqw = CQw[upper[1]];
                } else if (has_cq1 && has_cq0) {
                    Vector3d q_i = CQ.col(upper[0]);
                    Vector3d n_i = CQ.col(upper[0]);
                    Vector3d q_j = CQ.col(upper[1]);
                    Vector3d n_j = CQ.col(upper[1]);
                    auto result = mathext::compat_orientation_extrinsic_4(q_i, n_i, q_j, n_j);
                    cq = result.first * CQw[upper[0]] + result.second * CQw[upper[1]];
                    cqw = (CQw[upper[0]] + CQw[upper[1]]);
                }
                if (cq != Vector3d::Zero()) {
                    Vector3d n = N_next.col(i);
                    cq -= n.dot(cq) * n;
                    if (cq.squaredNorm() > RCPOVERFLOW) cq.normalize();
                }

                if (has_co0 && !has_co1) {
                    co = CO.col(upper[0]);
                    cow = COw[upper[0]];
                } else if (has_co1 && !has_co0) {
                    co = CO.col(upper[1]);
                    cow = COw[upper[1]];
                } else if (has_co1 && has_co0) {
                    double scale_x = m_scale;
                    double scale_y = m_scale;
                    if (with_scale) {
                        // FIXME
                        // scale_x *= S(0, i);
                        // scale_y *= S(1, i);
                    }
                    double inv_scale_x = 1.0f / scale_x;
                    double inv_scale_y = 1.0f / scale_y;

                    double scale_x_1 = m_scale;
                    double scale_y_1 = m_scale;
                    if (with_scale) {
                        // FIXME
                        // scale_x_1 *= S(0, j);
                        // scale_y_1 *= S(1, j);
                    }
                    double inv_scale_x_1 = 1.0f / scale_x_1;
                    double inv_scale_y_1 = 1.0f / scale_y_1;
                    auto result = mathext::compat_position_extrinsic_4(
                        V.col(upper[0]), N.col(upper[0]), CQ.col(upper[0]), CO.col(upper[0]),
                        V.col(upper[1]), N.col(upper[1]), CQ.col(upper[1]), CO.col(upper[1]), scale_x,
                        scale_y, inv_scale_x, inv_scale_y, scale_x_1, scale_y_1, inv_scale_x_1,
                        inv_scale_y_1);
                    cow = COw[upper[0]] + COw[upper[1]];
                    co = (result.first * COw[upper[0]] + result.second * COw[upper[1]]) / cow;
                }
                if (co != Vector3d::Zero()) {
                    Vector3d n = N_next.col(i), v = V_next.col(i);
                    co -= n.dot(cq - v) * n;
                }
                if (cqw > 0) cqw = 1;
                if (cow > 0) cow = 1;

                CQw_next[i] = cqw;
                COw_next[i] = cow;
                CQ_next.col(i) = cq;
                CO_next.col(i) = co;
            }
        }
    }
}


namespace quadriflow {
    using LinearSolver = Eigen::SimplicialLLT<Eigen::SparseMatrix<double> >;

    void Optimizer::optimize_positions(
        Hierarchy &hierarchy,
        const int iterations,
        bool with_scale
    ) {
        spdlog::info("Optimizing position field");

        const int n_adjacencies = static_cast<int>(hierarchy.m_normals.size());

        for (int level = n_adjacencies - 1; level >= 0; --level) {
            for (int iter = 0; iter < iterations; ++iter) {
                AdjacentMatrix &adj = hierarchy.m_adjacency[level];
                const MatrixXd &N = hierarchy.m_normals[level];
                const MatrixXd &Q = hierarchy.m_orientation[level];
                const MatrixXd &V = hierarchy.m_vertices[level];
                const MatrixXd &CQ = hierarchy.m_orientation_constraint[level];
                const MatrixXd &CO = hierarchy.m_position_constraints[level];
                const VectorXd &COw = hierarchy.m_position_constraint_weights[level];
                MatrixXd &O = hierarchy.m_positions[level];
                MatrixXd &S = hierarchy.m_scales[level];
                auto &phases = hierarchy.m_phases[level];
                for (int phase = 0; phase < phases.size(); ++phase) {
                    auto &p = phases[phase];
                    for (int pi = 0; pi < p.size(); ++pi) {
                        int i = p[pi];
                        double scale_x = hierarchy.m_scale;
                        double scale_y = hierarchy.m_scale;
                        if (with_scale) {
                            scale_x *= S(0, i);
                            scale_y *= S(1, i);
                        }
                        double inv_scale_x = 1.0f / scale_x;
                        double inv_scale_y = 1.0f / scale_y;
                        const Vector3d n_i = N.col(i), v_i = V.col(i);
                        Vector3d q_i = Q.col(i);

                        Vector3d sum = O.col(i);
                        double weight_sum = 0.0f;

                        q_i.normalize();
                        for (auto &link: adj[i]) {
                            const int j = link.id;
                            const double weight = link.weight;
                            if (weight == 0) continue;
                            double scale_x_1 = hierarchy.m_scale;
                            double scale_y_1 = hierarchy.m_scale;
                            if (with_scale) {
                                scale_x_1 *= S(0, j);
                                scale_y_1 *= S(1, j);
                            }
                            double inv_scale_x_1 = 1.0f / scale_x_1;
                            double inv_scale_y_1 = 1.0f / scale_y_1;

                            const Vector3d n_j = N.col(j), v_j = V.col(j);
                            Vector3d q_j = Q.col(j), o_j = O.col(j);

                            q_j.normalize();

                            std::pair<Vector3d, Vector3d> value = mathext::compat_position_extrinsic_4(
                                v_i, n_i, q_i, sum, v_j, n_j, q_j, o_j, scale_x, scale_y, inv_scale_x,
                                inv_scale_y, scale_x_1, scale_y_1, inv_scale_x_1, inv_scale_y_1);

                            sum = value.first * weight_sum + value.second * weight;
                            weight_sum += weight;
                            if (weight_sum > RCPOVERFLOW) sum /= weight_sum;
                            sum -= n_i.dot(sum - v_i) * n_i;
                        }

                        // Apply constraints
                        if (COw.size() > 0) {
                            float cw = COw[i];
                            if (cw != 0) {
                                Vector3d co = CO.col(i);
                                Vector3d cq = CQ.col(i);
                                Vector3d d = co - sum;
                                d -= cq.dot(d) * cq;
                                sum += cw * d;
                                sum -= n_i.dot(sum - v_i) * n_i;
                            }
                        }

                        if (weight_sum > 0) {
                            O.col(i) = mathext::position_round_4(sum, q_i, n_i, v_i, scale_x, scale_y,
                                                                 inv_scale_x, inv_scale_y);
                        }
                    }
                }
            }
            if (level > 0) {
                const MatrixXd &srcField = hierarchy.m_positions[level];
                const MatrixXi &toUpper = hierarchy.mToUpper[level - 1];
                MatrixXd &destField = hierarchy.m_positions[level - 1];
                const MatrixXd &N = hierarchy.m_normals[level - 1];
                const MatrixXd &V = hierarchy.m_vertices[level - 1];
                for (int i = 0; i < srcField.cols(); ++i) {
                    for (int k = 0; k < 2; ++k) {
                        int dest = toUpper(k, i);
                        if (dest == -1) continue;
                        Vector3d o = srcField.col(i), n = N.col(dest), v = V.col(dest);
                        o -= n * n.dot(o - v);
                        destField.col(dest) = o;
                    }
                }
            }
        }
    }

    void Optimizer::optimize_orientations(Hierarchy &hierarchy) {
        spdlog::info("Optimize orientations");

        const int levelIterations = 6;
        int n_normals = static_cast<int>(hierarchy.m_normals.size());
        for (int level = n_normals - 1; level >= 0; --level) {
            AdjacentMatrix &adj = hierarchy.m_adjacency[level];
            const MatrixXd &N = hierarchy.m_normals[level];
            const MatrixXd &CQ = hierarchy.m_orientation_constraint[level];
            const VectorXd &CQw = hierarchy.m_orientation_constraint_weight[level];
            MatrixXd &Q = hierarchy.m_orientation[level];
            auto &phases = hierarchy.m_phases[level];
            for (int iter = 0; iter < levelIterations; ++iter) {
                for (auto &p: phases) {
                    for (int i: p) {
                        const Vector3d n_i = N.col(i);
                        double weight_sum = 0.0f;
                        Vector3d sum = Q.col(i);
                        for (auto &link: adj[i]) {
                            const int j = link.id;
                            const double weight = link.weight;
                            if (weight == 0) continue;
                            const Vector3d n_j = N.col(j);
                            Vector3d q_j = Q.col(j);
                            std::pair<Vector3d, Vector3d> value =
                                    mathext::compat_orientation_extrinsic_4(sum, n_i, q_j, n_j);
                            sum = value.first * weight_sum + value.second * weight;
                            sum -= n_i * n_i.dot(sum);
                            weight_sum += weight;
                            double norm = sum.norm();
                            if (norm > RCPOVERFLOW) sum /= norm;
                        }

                        if (CQw.size() > 0) {
                            auto cw = CQw[i];
                            if (cw != 0) {
                                std::pair<Vector3d, Vector3d> value =
                                        mathext::compat_orientation_extrinsic_4(sum, n_i, CQ.col(i), n_i);
                                sum = value.first * (1 - cw) + value.second * cw;
                                sum -= n_i * n_i.dot(sum);

                                auto norm = sum.norm();
                                if (norm > RCPOVERFLOW) sum /= norm;
                            }
                        }

                        if (weight_sum > 0) {
                            Q.col(i) = sum;
                        }
                    }
                }
            }
            if (level > 0) {
                const MatrixXd &srcField = hierarchy.m_orientation[level];
                const MatrixXi &toUpper = hierarchy.mToUpper[level - 1];
                MatrixXd &destField = hierarchy.m_orientation[level - 1];
                const MatrixXd &normals_l1 = hierarchy.m_normals[level - 1];
                for (int i = 0; i < srcField.cols(); ++i) {
                    for (int k = 0; k < 2; ++k) {
                        int dest = toUpper(k, i);
                        if (dest == -1) continue;
                        Vector3d q = srcField.col(i), n = normals_l1.col(dest);
                        destField.col(dest) = q - n * n.dot(q);
                    }
                }
            }
        }

        for (int l = 0; l < hierarchy.m_normals.size() - 1; ++l) {
            const MatrixXd &N = hierarchy.m_normals[l];
            const MatrixXd &N_next = hierarchy.m_normals[l + 1];
            const MatrixXd &Q = hierarchy.m_orientation[l];
            MatrixXd &Q_next = hierarchy.m_orientation[l + 1];
            auto &toUpper = hierarchy.mToUpper[l];
            for (int i = 0; i < toUpper.cols(); ++i) {
                Vector2i upper = toUpper.col(i);
                Vector3d q0 = Q.col(upper[0]);
                Vector3d n0 = N.col(upper[0]);
                Vector3d q;

                if (upper[1] != -1) {
                    Vector3d q1 = Q.col(upper[1]);
                    Vector3d n1 = N.col(upper[1]);
                    auto result = mathext::compat_orientation_extrinsic_4(q0, n0, q1, n1);
                    q = result.first + result.second;
                } else {
                    q = q0;
                }
                Vector3d n = N_next.col(i);
                q -= n.dot(q) * n;
                if (q.squaredNorm() > RCPOVERFLOW) q.normalize();

                Q_next.col(i) = q;
            }
        }
    }

    void Optimizer::optimize_scale(Hierarchy &mRes, VectorXd &rho, int adaptive) {
        const MatrixXd &N = mRes.m_normals[0];
        MatrixXd &Q = mRes.m_orientation[0];
        MatrixXd &V = mRes.m_vertices[0];
        MatrixXd &S = mRes.m_scales[0];
        MatrixXd &K = mRes.m_areas[0];
        MatrixXi &F = mRes.m_faces;

        if (adaptive) {
            std::vector<Eigen::Triplet<double> > lhsTriplets;

            lhsTriplets.reserve(F.cols() * 6);
            for (int i = 0; i < V.cols(); ++i) {
                for (int j = 0; j < 2; ++j) {
                    S(j, i) = 1.0;
                    double sc1 = std::max(0.75 * S(j, i), rho[i] * 1.0 / mRes.m_scale);
                    S(j, i) = std::min(S(j, i), sc1);
                }
            }

            std::vector<std::map<int, double> > entries(V.cols() * 2);
            double lambda = 1;
            for (int i = 0; i < entries.size(); ++i) {
                entries[i][i] = lambda;
            }
            for (int i = 0; i < F.cols(); ++i) {
                for (int j = 0; j < 3; ++j) {
                    int v1 = F(j, i);
                    int v2 = F((j + 1) % 3, i);
                    Vector3d diff = V.col(v2) - V.col(v1);
                    Vector3d q_1 = Q.col(v1);
                    Vector3d q_2 = Q.col(v2);
                    Vector3d n_1 = N.col(v1);
                    Vector3d n_2 = N.col(v2);
                    Vector3d q_1_y = n_1.cross(q_1);
                    auto index = mathext::compat_orientation_extrinsic_index_4(q_1, n_1, q_2, n_2);
                    int v1_x = v1 * 2, v1_y = v1 * 2 + 1, v2_x = v2 * 2, v2_y = v2 * 2 + 1;

                    double dx = diff.dot(q_1);
                    double dy = diff.dot(q_1_y);

                    double kx_g = K(0, v1);
                    double ky_g = K(1, v1);

                    if (index.first % 2 != index.second % 2) {
                        std::swap(v2_x, v2_y);
                    }
                    double scale_x = (fmin(fmax(1 + kx_g * dy, 0.3), 3));
                    double scale_y = (fmin(fmax(1 + ky_g * dx, 0.3), 3));
                    //                (v2_x - scale_x * v1_x)^2 = 0
                    // x^2 - 2s xy + s^2 y^2
                    entries[v2_x][v2_x] += 1;
                    entries[v1_x][v1_x] += scale_x * scale_x;
                    entries[v2_y][v2_y] += 1;
                    entries[v1_y][v1_y] += scale_y * scale_y;
                    auto it = entries[v1_x].find(v2_x);
                    if (it == entries[v1_x].end()) {
                        entries[v1_x][v2_x] = -scale_x;
                        entries[v2_x][v1_x] = -scale_x;
                        entries[v1_y][v2_y] = -scale_y;
                        entries[v2_y][v1_y] = -scale_y;
                    } else {
                        it->second -= scale_x;
                        entries[v2_x][v1_x] -= scale_x;
                        entries[v1_y][v2_y] -= scale_y;
                        entries[v2_y][v1_y] -= scale_y;
                    }
                }
            }

            Eigen::SparseMatrix<double> A(V.cols() * 2, V.cols() * 2);
            VectorXd rhs(V.cols() * 2);
            rhs.setZero();
            for (int i = 0; i < entries.size(); ++i) {
                rhs(i) = lambda * S(i % 2, i / 2);
                for (auto &rec: entries[i]) {
                    lhsTriplets.push_back(Eigen::Triplet<double>(i, rec.first, rec.second));
                }
            }
            A.setFromTriplets(lhsTriplets.begin(), lhsTriplets.end());
            LinearSolver solver;
            solver.analyzePattern(A);

            solver.factorize(A);

            VectorXd result = solver.solve(rhs);

            double total_area = 0;
            for (int i = 0; i < V.cols(); ++i) {
                S(0, i) = (result(i * 2));
                S(1, i) = (result(i * 2 + 1));
                total_area += S(0, i) * S(1, i);
            }
            total_area = sqrt(V.cols() / total_area);
            for (int i = 0; i < V.cols(); ++i) {
                //            S(0, i) *= total_area;
                //            S(1, i) *= total_area;
            }
        } else {
            for (int i = 0; i < V.cols(); ++i) {
                S(0, i) = 1;
                S(1, i) = 1;
            }
        }

        for (int l = 0; l < mRes.m_scales.size() - 1; ++l) {
            const MatrixXd &S = mRes.m_scales[l];
            MatrixXd &S_next = mRes.m_scales[l + 1];
            auto &toUpper = mRes.mToUpper[l];
            for (int i = 0; i < toUpper.cols(); ++i) {
                Vector2i upper = toUpper.col(i);
                Vector2d q0 = S.col(upper[0]);

                if (upper[1] != -1) {
                    q0 = (q0 + S.col(upper[1])) * 0.5;
                }
                S_next.col(i) = q0;
            }
        }
    }

    void Optimizer::optimize_positions_dynamic(
        MatrixXi &F,
        MatrixXd &V,
        MatrixXd &N,
        MatrixXd &Q,
        std::vector<std::vector<int> > &Vset,
        std::vector<Vector3d> &O_compact,
        std::vector<Vector4i> &F_compact,
        std::vector<int> &V2E_compact,
        std::vector<int> &E2E_compact,
        double mScale,
        std::vector<Vector3d> &diffs,
        std::vector<int> &diff_count,
        std::map<std::pair<int, int>, int> &o2e,
        std::vector<int> &sharp_o,
        std::map<int, std::pair<Vector3d, Vector3d> > &compact_sharp_constraints,
        bool with_scale
    ) {
        std::set<int> uncertain;
        for (auto &info: o2e) {
            if (diff_count[info.second] == 0) {
                uncertain.insert(info.first.first);
                uncertain.insert(info.first.second);
            }
        }
        std::vector<int> Vind(O_compact.size(), -1);
        std::vector<std::list<int> > links(O_compact.size());
        std::vector<std::list<int> > dedges(O_compact.size());
        std::vector<std::vector<int> > adj(V.cols());
        for (int i = 0; i < F.cols(); ++i) {
            for (int j = 0; j < 3; ++j) {
                int v1 = F(j, i);
                int v2 = F((j + 1) % 3, i);
                adj[v1].push_back(v2);
            }
        }
        auto FindNearest = [&]() {
            for (int i = 0; i < O_compact.size(); ++i) {
                if (Vind[i] == -1) {
                    double min_dis = 1e30;
                    int min_ind = -1;
                    for (auto v: Vset[i]) {
                        double dis = (V.col(v) - O_compact[i]).squaredNorm();
                        if (dis < min_dis) {
                            min_dis = dis;
                            min_ind = v;
                        }
                    }
                    if (min_ind > -1) {
                        Vind[i] = min_ind;
                        double x = (O_compact[i] - V.col(min_ind)).dot(N.col(min_ind));
                        O_compact[i] -= x * N.col(min_ind);
                    }
                } else {
                    int current_v = Vind[i];
                    Vector3d n = N.col(current_v);
                    double current_dis = (O_compact[i] - V.col(current_v)).squaredNorm();
                    while (true) {
                        int next_v = -1;
                        for (auto &v: adj[current_v]) {
                            if (N.col(v).dot(n) < cos(10.0 / 180.0 * 3.141592654)) continue;
                            double dis = (O_compact[i] - V.col(v)).squaredNorm();
                            if (dis < current_dis) {
                                current_dis = dis;
                                next_v = v;
                            }
                        }
                        if (next_v == -1) break;
                        // rotate ideal distance
                        Vector3d n1 = N.col(current_v);
                        Vector3d n2 = N.col(next_v);
                        Vector3d axis = n1.cross(n2);
                        double len = axis.norm();
                        double angle = atan2(len, n1.dot(n2));
                        axis.normalized();
                        Matrix3d m = AngleAxisd(angle, axis).toRotationMatrix();
                        for (auto e: dedges[i]) {
                            Vector3d &d = diffs[e];
                            d = m * d;
                        }
                        current_v = next_v;
                    }
                    Vind[i] = current_v;
                }
            }
        };

        auto BuildConnection = [&]() {
            for (int i = 0; i < links.size(); ++i) {
                int deid0 = V2E_compact[i];
                if (deid0 != -1) {
                    std::list<int> &connection = links[i];
                    std::list<int> &dedge = dedges[i];
                    int deid = deid0;
                    do {
                        connection.push_back(F_compact[deid / 4][(deid + 1) % 4]);
                        dedge.push_back(deid);
                        deid = E2E_compact[deid / 4 * 4 + (deid + 3) % 4];
                    } while (deid != -1 && deid != deid0);
                    if (deid == -1) {
                        deid = deid0;
                        do {
                            deid = E2E_compact[deid];
                            if (deid == -1) break;
                            deid = deid / 4 * 4 + (deid + 1) % 4;
                            connection.push_front(F_compact[deid / 4][(deid + 1) % 4]);
                            dedge.push_front(deid);
                        } while (true);
                    }
                }
            }
        };

        std::vector<Vector3d> lines;
        auto ComputeDistance = [&]() {
            std::set<int> unobserved;
            for (auto &info: o2e) {
                if (diff_count[info.second] == 0) {
                    unobserved.insert(info.first.first);
                }
            }
            while (true) {
                bool update = false;
                std::set<int> observed;
                for (auto &p: unobserved) {
                    std::vector<int> observations, edges;
                    int count = 0;
                    for (auto &e: dedges[p]) {
                        edges.push_back(e);
                        if (diff_count[e]) {
                            count += 1;
                            observations.push_back(1);
                        } else {
                            observations.push_back(0);
                        }
                    }
                    if (count <= 1) continue;
                    update = true;
                    observed.insert(p);
                    for (int i = 0; i < observations.size(); ++i) {
                        if (observations[i] == 1) continue;
                        int j = i;
                        std::list<int> interp;
                        while (observations[j] == 0) {
                            interp.push_front(j);
                            j -= 1;
                            if (j < 0) j = edges.size() - 1;
                        }
                        j = (i + 1) % edges.size();
                        while (observations[j] == 0) {
                            interp.push_back(j);
                            j += 1;
                            if (j == edges.size()) j = 0;
                        }
                        Vector3d dl = diffs[edges[(interp.front() + edges.size() - 1) % edges.size()]];
                        double lenl = dl.norm();
                        Vector3d dr = diffs[edges[(interp.back() + 1) % edges.size()]];
                        double lenr = dr.norm();
                        dl /= lenl;
                        dr /= lenr;
                        Vector3d n = dl.cross(dr).normalized();
                        double angle = atan2(dl.cross(dr).norm(), dl.dot(dr));
                        if (angle < 0) angle += 2 * 3.141592654;
                        Vector3d nc = N.col(Vind[p]);
                        if (n.dot(nc) < 0) {
                            n = -n;
                            angle = 2 * 3.141592654 - angle;
                        }
                        double step = (lenr - lenl) / (interp.size() + 1);
                        angle /= interp.size() + 1;
                        Vector3d dlp = nc.cross(dl).normalized();
                        int t = 0;
                        for (auto q: interp) {
                            t += 1;
                            observations[q] = 1;
                            double ad = angle * t;
                            int e = edges[q];
                            int re = E2E_compact[e];
                            diff_count[e] = 2;
                            diffs[e] = (cos(ad) * dl + sin(ad) * dlp) * (lenl + step * t);
                            if (re != -1) {
                                diff_count[re] = 2;
                                diffs[re] = -diffs[e];
                            }
                        }
                        for (int i = 0; i < edges.size(); ++i) {
                            lines.push_back(O_compact[p]);
                            lines.push_back(O_compact[p] + diffs[edges[i]]);
                        }
                    }
                }
                if (!update) break;
                for (auto &p: observed) unobserved.erase(p);
            }
        };

        BuildConnection();
        int max_iter = 10;
        for (int iter = 0; iter < max_iter; ++iter) {
            FindNearest();
            ComputeDistance();

            std::vector<std::unordered_map<int, double> > entries(O_compact.size() * 2);
            std::vector<int> fixed_dim(O_compact.size() * 2, 0);
            for (auto &info: compact_sharp_constraints) {
                fixed_dim[info.first * 2 + 1] = 1;
                if (info.second.second.norm() < 0.5) fixed_dim[info.first * 2] = 1;
            }
            std::vector<double> b(O_compact.size() * 2);
            std::vector<double> x(O_compact.size() * 2);
            std::vector<Vector3d> Q_compact(O_compact.size());
            std::vector<Vector3d> N_compact(O_compact.size());
            std::vector<Vector3d> V_compact(O_compact.size());
            for (int i = 0; i < O_compact.size(); ++i) {
                Q_compact[i] = Q.col(Vind[i]);
                N_compact[i] = N.col(Vind[i]);
                V_compact[i] = V.col(Vind[i]);
                if (fixed_dim[i * 2 + 1] && !fixed_dim[i * 2]) {
                    Q_compact[i] = compact_sharp_constraints[i].second;
                    V_compact[i] = compact_sharp_constraints[i].first;
                }
            }
            for (int i = 0; i < O_compact.size(); ++i) {
                Vector3d q = Q_compact[i];
                Vector3d n = N_compact[i];
                Vector3d q_y = n.cross(q);
                auto Vi = V_compact[i];
                x[i * 2] = (O_compact[i] - Vi).dot(q);
                x[i * 2 + 1] = (O_compact[i] - Vi).dot(q_y);
            }
            for (int i = 0; i < O_compact.size(); ++i) {
                Vector3d qx = Q_compact[i];
                Vector3d qy = N_compact[i];
                qy = qy.cross(qx);
                auto dedge_it = dedges[i].begin();
                for (auto it = links[i].begin(); it != links[i].end(); ++it, ++dedge_it) {
                    int j = *it;
                    Vector3d qx2 = Q_compact[j];
                    Vector3d qy2 = N_compact[j];
                    qy2 = qy2.cross(qx2);

                    int de = o2e[std::make_pair(i, j)];
                    double lambda = (diff_count[de] == 1) ? 1 : 1;
                    Vector3d target_offset = diffs[de];

                    auto Vi = V_compact[i];
                    auto Vj = V_compact[j];

                    Vector3d offset = Vj - Vi;

                    //                target_offset.normalize();
                    //                target_offset *= m_scale;
                    Vector3d C = target_offset - offset;
                    int vid[] = {j * 2, j * 2 + 1, i * 2, i * 2 + 1};
                    Vector3d weights[] = {qx2, qy2, -qx, -qy};
                    for (int ii = 0; ii < 4; ++ii) {
                        for (int jj = 0; jj < 4; ++jj) {
                            auto it = entries[vid[ii]].find(vid[jj]);
                            if (it == entries[vid[ii]].end()) {
                                entries[vid[ii]][vid[jj]] = lambda * weights[ii].dot(weights[jj]);
                            } else {
                                entries[vid[ii]][vid[jj]] += lambda * weights[ii].dot(weights[jj]);
                            }
                        }
                        b[vid[ii]] += lambda * weights[ii].dot(C);
                    }
                }
            }

            // fix sharp edges
            for (int i = 0; i < entries.size(); ++i) {
                if (entries[i].size() == 0) {
                    entries[i][i] = 1;
                    b[i] = x[i];
                }
                if (fixed_dim[i]) {
                    b[i] = x[i];
                    entries[i].clear();
                    entries[i][i] = 1;
                } else {
                    std::unordered_map<int, double> newmap;
                    for (auto &rec: entries[i]) {
                        if (fixed_dim[rec.first]) {
                            b[i] -= rec.second * x[rec.first];
                        } else {
                            newmap[rec.first] = rec.second;
                        }
                    }
                    std::swap(entries[i], newmap);
                }
            }
            std::vector<Eigen::Triplet<double> > lhsTriplets;
            lhsTriplets.reserve(F_compact.size() * 8);
            Eigen::SparseMatrix<double> A(O_compact.size() * 2, O_compact.size() * 2);
            VectorXd rhs(O_compact.size() * 2);
            rhs.setZero();
            for (int i = 0; i < entries.size(); ++i) {
                rhs(i) = b[i];
                for (auto &rec: entries[i]) {
                    lhsTriplets.push_back(Eigen::Triplet<double>(i, rec.first, rec.second));
                }
            }

            A.setFromTriplets(lhsTriplets.begin(), lhsTriplets.end());

            // FIXME: IncompleteCholesky Preconditioner will fail here so I fallback to Diagonal one.
            // I suspected either there is a implementation bug in IncompleteCholesky Preconditioner
            // or there is a memory corruption somewhere.  However, g++'s address sanitizer does not
            // report anything useful.
            LinearSolver solver;
            solver.analyzePattern(A);
            solver.factorize(A);
            //        Eigen::setNbThreads(1);
            //        ConjugateGradient<SparseMatrix<double>, Lower | Upper> solver;
            //        VectorXd x0 = VectorXd::Map(x.data(), x.size());
            //        solver.setMaxIterations(40);

            //        solver.compute(m_vertex_area);
            VectorXd x_new = solver.solve(rhs); // solver.solveWithGuess(rhs, x0);

            //            int t2 = GetCurrentTime64();
            //            printf("[LSQ] Linear solver uses %lf seconds.\n", (t2 - t1) * 1e-3);
            for (int i = 0; i < O_compact.size(); ++i) {
                // Vector3d q = Q.col(Vind[i]);
                Vector3d q = Q_compact[i];
                // Vector3d n = m_normals_vertices.col(Vind[i]);
                Vector3d n = N_compact[i];
                Vector3d q_y = n.cross(q);
                auto Vi = V_compact[i];
                O_compact[i] = Vi + q * x_new[i * 2] + q_y * x_new[i * 2 + 1];
            }

            // forgive my hack...
            if (iter + 1 == max_iter) {
                for (int iter = 0; iter < 5; ++iter) {
                    for (int i = 0; i < O_compact.size(); ++i) {
                        if (sharp_o[i]) continue;
                        if (dedges[i].size() != 4 || uncertain.count(i)) {
                            Vector3d n(0, 0, 0), v(0, 0, 0);
                            Vector3d v0 = O_compact[i];
                            for (auto e: dedges[i]) {
                                Vector3d v1 = O_compact[F_compact[e / 4][(e + 1) % 4]];
                                Vector3d v2 = O_compact[F_compact[e / 4][(e + 3) % 4]];
                                n += (v1 - v0).cross(v2 - v0);
                                v += v1;
                            }
                            n.normalize();
                            Vector3d offset = v / dedges[i].size() - v0;
                            offset -= offset.dot(n) * n;
                            O_compact[i] += offset;
                        }
                    }
                }
            }
        }
    }

    void Optimizer::optimize_positions_sharp(
        Hierarchy &mRes,
        std::vector<DEdge> &edge_values,
        std::vector<Vector2i> &edge_diff,
        std::vector<int> &sharp_edges,
        std::set<int> &sharp_vertices,
        std::map<int, std::pair<Vector3d, Vector3d> > &sharp_constraints,
        bool with_scale
    ) {
        spdlog::debug("Optimizing edge positions sharp_edges={}, sharp_vertices={}", sharp_edges.size(), sharp_vertices.size());

        auto &V = mRes.m_vertices[0];
        auto &F = mRes.m_faces;
        auto &Q = mRes.m_orientation[0];
        auto &N = mRes.m_normals[0];
        auto &O = mRes.m_positions[0];
        auto &S = mRes.m_scales[0];

        DisjointTree tree(V.cols());
        for (int i = 0; i < edge_diff.size(); ++i) {
            if (edge_diff[i].array().abs().sum() == 0) {
                tree.Merge(edge_values[i].x, edge_values[i].y);
            }
        }
        tree.BuildCompactParent();
        std::map<int, int> compact_sharp_indices;
        std::set<DEdge> compact_sharp_edges;
        for (int i = 0; i < sharp_edges.size(); ++i) {
            if (sharp_edges[i] == 1) {
                int v1 = tree.Index(F(i % 3, i / 3));
                int v2 = tree.Index(F((i + 1) % 3, i / 3));
                compact_sharp_edges.insert(DEdge(v1, v2));
            }
        }
        for (auto &v: sharp_vertices) {
            int p = tree.Index(v);
            if (compact_sharp_indices.count(p) == 0) {
                int s = compact_sharp_indices.size();
                compact_sharp_indices[p] = s;
            }
        }
        std::map<int, std::set<int> > sharp_vertices_links;
        std::set<DEdge> sharp_dedges;
        for (int i = 0; i < sharp_edges.size(); ++i) {
            if (sharp_edges[i]) {
                int v1 = F(i % 3, i / 3);
                int v2 = F((i + 1) % 3, i / 3);
                if (sharp_vertices_links.count(v1) == 0) sharp_vertices_links[v1] = std::set<int>();
                sharp_vertices_links[v1].insert(v2);
                sharp_dedges.insert(DEdge(v1, v2));
            }
        }
        std::vector<std::vector<int> > sharp_to_original_indices(compact_sharp_indices.size());
        for (auto &v: sharp_vertices_links) {
            if (v.second.size() == 2) continue;
            int p = tree.Index(v.first);
            sharp_to_original_indices[compact_sharp_indices[p]].push_back(v.first);
        }
        for (auto &v: sharp_vertices_links) {
            if (v.second.size() != 2) continue;
            int p = tree.Index(v.first);
            sharp_to_original_indices[compact_sharp_indices[p]].push_back(v.first);
        }

        for (int i = 0; i < V.cols(); ++i) {
            if (sharp_vertices.count(i)) continue;
            int p = tree.Index(i);
            if (compact_sharp_indices.count(p))
                sharp_to_original_indices[compact_sharp_indices[p]].push_back(i);
        }

        int num = sharp_to_original_indices.size();
        std::vector<std::set<int> > links(sharp_to_original_indices.size());
        for (int e = 0; e < edge_diff.size(); ++e) {
            int v1 = edge_values[e].x;
            int v2 = edge_values[e].y;
            int p1 = tree.Index(v1);
            int p2 = tree.Index(v2);
            if (p1 == p2 || compact_sharp_edges.count(DEdge(p1, p2)) == 0) continue;
            p1 = compact_sharp_indices[p1];
            p2 = compact_sharp_indices[p2];

            links[p1].insert(p2);
            links[p2].insert(p1);
        }

        std::vector<int> hash(links.size(), 0);
        std::vector<std::vector<Vector3d> > loops;
        for (int i = 0; i < num; ++i) {
            if (hash[i] == 1) continue;
            if (links[i].size() == 2) {
                std::vector<int> q;
                q.push_back(i);
                hash[i] = 1;
                int v = i;
                int prev_v = -1;
                bool is_loop = false;
                while (links[v].size() == 2) {
                    int next_v = -1;
                    for (auto nv: links[v])
                        if (nv != prev_v) next_v = nv;
                    if (hash[next_v]) {
                        is_loop = true;
                        break;
                    }
                    if (links[next_v].size() == 2) hash[next_v] = true;
                    q.push_back(next_v);
                    prev_v = v;
                    v = next_v;
                }
                if (!is_loop && q.size() >= 2) {
                    std::vector<int> q1;
                    int v = i;
                    int prev_v = q[1];
                    while (links[v].size() == 2) {
                        int next_v = -1;
                        for (auto nv: links[v])
                            if (nv != prev_v) next_v = nv;
                        if (hash[next_v]) {
                            is_loop = true;
                            break;
                        }
                        if (links[next_v].size() == 2) hash[next_v] = true;
                        q1.push_back(next_v);
                        prev_v = v;
                        v = next_v;
                    }
                    std::reverse(q1.begin(), q1.end());
                    q1.insert(q1.end(), q.begin(), q.end());
                    std::swap(q1, q);
                }
                if (q.size() < 3) continue;
                if (is_loop) q.push_back(q.front());
                double len = 0, scale = 0;
                std::vector<Vector3d> o(q.size()), new_o(q.size());
                std::vector<double> sc(q.size());

                for (int i = 0; i < q.size() - 1; ++i) {
                    int v1 = q[i];
                    int v2 = q[i + 1];
                    auto it = links[v1].find(v2);
                    if (it == links[v1].end()) {
                        printf("Non exist!\n");
                        exit(0);
                    }
                }

                for (int i = 0; i < q.size(); ++i) {
                    if (sharp_to_original_indices[q[i]].size() == 0) {
                        continue;
                    }
                    o[i] = O.col(sharp_to_original_indices[q[i]][0]);
                    Vector3d qx = Q.col(sharp_to_original_indices[q[i]][0]);
                    Vector3d qy = Vector3d(N.col(sharp_to_original_indices[q[i]][0])).cross(qx);
                    int fst = sharp_to_original_indices[q[1]][0];
                    Vector3d dis = (i == 0) ? (Vector3d(O.col(fst)) - o[i]) : o[i] - o[i - 1];
                    if (with_scale)
                        sc[i] = (abs(qx.dot(dis)) > abs(qy.dot(dis)))
                                    ? S(0, sharp_to_original_indices[q[i]][0])
                                    : S(1, sharp_to_original_indices[q[i]][0]);
                    else
                        sc[i] = 1;
                    new_o[i] = o[i];
                }

                if (is_loop) {
                    for (int i = 0; i < q.size(); ++i) {
                        Vector3d dir =
                                (o[(i + 1) % q.size()] - o[(i + q.size() - 1) % q.size()]).normalized();
                        for (auto &ind: sharp_to_original_indices[q[i]]) {
                            sharp_constraints[ind] = std::make_pair(o[i], dir);
                        }
                    }
                } else {
                    for (int i = 0; i < q.size(); ++i) {
                        Vector3d dir(0, 0, 0);
                        if (i != 0 && i + 1 != q.size())
                            dir = (o[i + 1] - o[i - 1]).normalized();
                        else if (links[q[i]].size() == 1) {
                            if (i == 0)
                                dir = (o[i + 1] - o[i]).normalized();
                            else
                                dir = (o[i] - o[i - 1]).normalized();
                        }
                        for (auto &ind: sharp_to_original_indices[q[i]]) {
                            sharp_constraints[ind] = std::make_pair(o[i], dir);
                        }
                    }
                }

                for (int i = 0; i < q.size() - 1; ++i) {
                    len += (o[i + 1] - o[i]).norm();
                    scale += sc[i];
                }

                int next_m = q.size() - 1;

                double left_norm = len * sc[0] / scale;
                int current_v = 0;
                double current_norm = (o[1] - o[0]).norm();
                for (int i = 1; i < next_m; ++i) {
                    while (left_norm >= current_norm) {
                        left_norm -= current_norm;
                        current_v += 1;
                        current_norm = (o[current_v + 1] - o[current_v]).norm();
                    }
                    new_o[i] =
                            (o[current_v + 1] * left_norm + o[current_v] * (current_norm - left_norm)) /
                            current_norm;
                    o[current_v] = new_o[i];
                    current_norm -= left_norm;
                    left_norm = len * sc[current_v] / scale;
                }

                for (int i = 0; i < q.size(); ++i) {
                    for (auto v: sharp_to_original_indices[q[i]]) {
                        O.col(v) = new_o[i];
                    }
                }

                loops.push_back(new_o);
            }
        }
    }

    void Optimizer::optimize_positions_fixed(
        Hierarchy &mRes,
        std::vector<DEdge> &edge_values,
        std::vector<Vector2i> &edge_diff,
        std::set<int> &sharp_vertices,
        std::map<int, std::pair<Vector3d, Vector3d> > &sharp_constraints,
        bool with_scale
    ) {
        spdlog::debug("Optimizing fixed vertex positions sharp_vertices={}", sharp_vertices.size());

        auto &V = mRes.m_vertices[0];
        auto &F = mRes.m_faces;
        auto &Q = mRes.m_orientation[0];
        auto &N = mRes.m_normals[0];
        auto &O = mRes.m_positions[0];
        auto &S = mRes.m_scales[0];

        DisjointTree tree(V.cols());
        for (int i = 0; i < edge_diff.size(); ++i) {
            if (edge_diff[i].array().abs().sum() == 0) {
                tree.Merge(edge_values[i].x, edge_values[i].y);
            }
        }
        tree.BuildCompactParent();
        int num = tree.CompactNum();

        // Find the most descriptive vertex
        std::vector<Vector3d> v_positions(num, Vector3d(0, 0, 0));
        std::vector<int> v_count(num);
        std::vector<double> v_distance(num, 1e30);
        std::vector<int> v_index(num, -1);

        for (int i = 0; i < V.cols(); ++i) {
            v_positions[tree.Index(i)] += O.col(i);
            v_count[tree.Index(i)] += 1;
        }
        for (int i = 0; i < num; ++i) {
            if (v_count[i] > 0) v_positions[i] /= v_count[i];
        }
        for (int i = 0; i < V.cols(); ++i) {
            int p = tree.Index(i);
            double dis = (v_positions[p] - V.col(i)).squaredNorm();
            if (dis < v_distance[p]) {
                v_distance[p] = dis;
                v_index[p] = i;
            }
        }

        std::set<int> compact_sharp_vertices;
        for (auto &v: sharp_vertices) {
            v_positions[tree.Index(v)] = O.col(v);
            v_index[tree.Index(v)] = v;
            V.col(v) = O.col(v);
            compact_sharp_vertices.insert(tree.Index(v));
        }
        std::vector<std::map<int, std::pair<int, Vector3d> > > ideal_distances(tree.CompactNum());
        for (int e = 0; e < edge_diff.size(); ++e) {
            int v1 = edge_values[e].x;
            int v2 = edge_values[e].y;

            int p1 = tree.Index(v1);
            int p2 = tree.Index(v2);
            int q1 = v_index[p1];
            int q2 = v_index[p2];

            Vector3d q_1 = Q.col(v1);
            Vector3d q_2 = Q.col(v2);

            Vector3d n_1 = N.col(v1);
            Vector3d n_2 = N.col(v2);
            Vector3d q_1_y = n_1.cross(q_1);
            Vector3d q_2_y = n_2.cross(q_2);
            auto index = mathext::compat_orientation_extrinsic_index_4(q_1, n_1, q_2, n_2);
            double s_x1 = S(0, v1), s_y1 = S(1, v1);
            double s_x2 = S(0, v2), s_y2 = S(1, v2);
            int rank_diff = (index.second + 4 - index.first) % 4;
            if (rank_diff % 2 == 1) std::swap(s_x2, s_y2);
            Vector3d qd_x = 0.5 * (mathext::rotate90_by(q_2, n_2, rank_diff) + q_1);
            Vector3d qd_y = 0.5 * (mathext::rotate90_by(q_2_y, n_2, rank_diff) + q_1_y);
            double scale_x = (with_scale ? 0.5 * (s_x1 + s_x2) : 1) * mRes.m_scale;
            double scale_y = (with_scale ? 0.5 * (s_y1 + s_y2) : 1) * mRes.m_scale;
            Vector2i diff = edge_diff[e];

            Vector3d origin1 =
                    /*(sharp_constraints.count(q1)) ? sharp_constraints[q1].first : */ V.col(q1);
            Vector3d origin2 =
                    /*(sharp_constraints.count(q2)) ? sharp_constraints[q2].first : */ V.col(q2);
            Vector3d C = diff[0] * scale_x * qd_x + diff[1] * scale_y * qd_y + origin1 - origin2;
            auto it = ideal_distances[p1].find(p2);
            if (it == ideal_distances[p1].end()) {
                ideal_distances[p1][p2] = std::make_pair(1, C);
            } else {
                it->second.first += 1;
                it->second.second += C;
            }
        }

        std::vector<std::unordered_map<int, double> > entries(num * 2);
        std::vector<double> b(num * 2);

        for (int m = 0; m < num; ++m) {
            int v1 = v_index[m];
            for (auto &info: ideal_distances[m]) {
                int v2 = v_index[info.first];
                Vector3d q_1 = Q.col(v1);
                Vector3d q_2 = Q.col(v2);
                if (sharp_constraints.count(v1)) {
                    Vector3d d = sharp_constraints[v1].second;
                    if (d != Vector3d::Zero()) q_1 = d;
                }
                if (sharp_constraints.count(v2)) {
                    Vector3d d = sharp_constraints[v2].second;
                    if (d != Vector3d::Zero()) q_2 = d;
                }

                Vector3d n_1 = N.col(v1);
                Vector3d n_2 = N.col(v2);
                Vector3d q_1_y = n_1.cross(q_1);
                Vector3d q_2_y = n_2.cross(q_2);
                Vector3d weights[] = {q_2, q_2_y, -q_1, -q_1_y};
                int vid[] = {info.first * 2, info.first * 2 + 1, m * 2, m * 2 + 1};
                Vector3d dis = info.second.second / info.second.first;
                double lambda = 1;
                if (sharp_vertices.count(v1) && sharp_vertices.count(v2)) lambda = 1;
                for (int i = 0; i < 4; ++i) {
                    for (int j = 0; j < 4; ++j) {
                        auto it = entries[vid[i]].find(vid[j]);
                        if (it == entries[vid[i]].end()) {
                            entries[vid[i]][vid[j]] = weights[i].dot(weights[j]) * lambda;
                        } else {
                            entries[vid[i]][vid[j]] += weights[i].dot(weights[j]) * lambda;
                        }
                    }
                    b[vid[i]] += weights[i].dot(dis) * lambda;
                }
            }
        }

        std::vector<int> fixed_dim(num * 2, 0);
        std::vector<double> x(num * 2);
        for (int i = 0; i < num; ++i) {
            int p = v_index[i];
            Vector3d q = Q.col(p);

            if (sharp_constraints.count(p)) {
                Vector3d dir = sharp_constraints[p].second;
                fixed_dim[i * 2 + 1] = 1;
                if (dir != Vector3d::Zero()) {
                    q = dir;
                } else
                    fixed_dim[i * 2] = 1;
            }
            Vector3d n = N.col(p);
            Vector3d q_y = n.cross(q);
            x[i * 2] = (v_positions[i] - V.col(p)).dot(q);
            x[i * 2 + 1] = (v_positions[i] - V.col(p)).dot(q_y);
        }

        // fix sharp edges
        for (int i = 0; i < entries.size(); ++i) {
            if (fixed_dim[i]) {
                b[i] = x[i];
                entries[i].clear();
                entries[i][i] = 1;
            } else {
                std::unordered_map<int, double> newmap;
                for (auto &rec: entries[i]) {
                    if (fixed_dim[rec.first]) {
                        b[i] -= rec.second * x[rec.first];
                    } else {
                        newmap[rec.first] = rec.second;
                    }
                }
                std::swap(entries[i], newmap);
            }
        }
        for (int i = 0; i < entries.size(); ++i) {
            if (entries[i].size() == 0) {
                entries[i][i] = 1;
            }
        }

        std::vector<Eigen::Triplet<double> > lhsTriplets;
        lhsTriplets.reserve(F.cols() * 6);
        Eigen::SparseMatrix<double> A(num * 2, num * 2);
        VectorXd rhs(num * 2);
        rhs.setZero();
        for (int i = 0; i < entries.size(); ++i) {
            rhs(i) = b[i];
            if (std::isnan(b[i])) {
                printf("Equation has nan!\n");
                exit(0);
            }
            for (auto &rec: entries[i]) {
                lhsTriplets.push_back(Eigen::Triplet<double>(i, rec.first, rec.second));
                if (std::isnan(rec.second)) {
                    printf("Equation has nan!\n");
                    exit(0);
                }
            }
        }
        A.setFromTriplets(lhsTriplets.begin(), lhsTriplets.end());

        LinearSolver solver;
        solver.analyzePattern(A);
        solver.factorize(A);

        VectorXd x_new = solver.solve(rhs);

        for (int i = 0; i < x.size(); ++i) {
            if (!std::isnan(x_new[i])) {
                if (!fixed_dim[i / 2 * 2 + 1]) {
                    double total = 0;
                    for (auto &res: entries[i]) {
                        double t = x_new[res.first];
                        if (std::isnan(t)) t = 0;
                        total += t * res.second;
                    }
                }
                x[i] = x_new[i];
            }
        }

        for (int i = 0; i < O.cols(); ++i) {
            int p = tree.Index(i);
            int c = v_index[p];
            Vector3d q = Q.col(c);
            if (fixed_dim[p * 2 + 1]) {
                Vector3d dir = sharp_constraints[c].second;
                if (dir != Vector3d::Zero()) q = dir;
            }
            Vector3d n = N.col(c);
            Vector3d q_y = n.cross(q);
            O.col(i) = V.col(c) + q * x[p * 2] + q_y * x[p * 2 + 1];
        }
    }

    void Optimizer::optimize_integer_constraints(Hierarchy &mRes, std::map<int, int> &singularities) {
        int edge_capacity = 2;
        bool fullFlow = false;
        std::vector<std::vector<int> > &AllowChange = mRes.mAllowChanges;
        for (int level = mRes.mToUpperEdges.size(); level >= 0; --level) {
            auto &EdgeDiff = mRes.mEdgeDiff[level];
            auto &FQ = mRes.mFQ[level];
            auto &F2E = mRes.mF2E[level];
            auto &E2F = mRes.mE2F[level];

            int iter = 0;
            while (!fullFlow) {
                std::vector<Vector4i> edge_to_constraints(E2F.size() * 2, Vector4i(-1, 0, -1, 0));
                std::vector<int> initial(F2E.size() * 2, 0);
                for (int i = 0; i < F2E.size(); ++i) {
                    for (int j = 0; j < 3; ++j) {
                        int e = F2E[i][j];
                        Vector2i index = mathext::rshift90(Vector2i(e * 2 + 1, e * 2 + 2), FQ[i][j]);
                        for (int k = 0; k < 2; ++k) {
                            int l = abs(index[k]);
                            int s = index[k] / l;
                            int ind = l - 1;
                            int equationID = i * 2 + k;
                            if (edge_to_constraints[ind][0] == -1) {
                                edge_to_constraints[ind][0] = equationID;
                                edge_to_constraints[ind][1] = s;
                            } else {
                                edge_to_constraints[ind][2] = equationID;
                                edge_to_constraints[ind][3] = s;
                            }
                            initial[equationID] += s * EdgeDiff[ind / 2][ind % 2];
                        }
                    }
                }
                std::vector<std::pair<Vector2i, int> > arcs;
                std::vector<int> arc_ids;
                for (int i = 0; i < edge_to_constraints.size(); ++i) {
                    if (AllowChange[level][i] == 0) continue;
                    if (edge_to_constraints[i][0] == -1 || edge_to_constraints[i][2] == -1) continue;
                    if (edge_to_constraints[i][1] == -edge_to_constraints[i][3]) {
                        int v1 = edge_to_constraints[i][0];
                        int v2 = edge_to_constraints[i][2];
                        if (edge_to_constraints[i][1] < 0) std::swap(v1, v2);
                        int current_v = EdgeDiff[i / 2][i % 2];
                        arcs.push_back(std::make_pair(Vector2i(v1, v2), current_v));
                        if (AllowChange[level][i] == 1)
                            arc_ids.push_back(i + 1);
                        else {
                            arc_ids.push_back(-(i + 1));
                        }
                    }
                }
                int supply = 0;
                int demand = 0;
                for (int i = 0; i < initial.size(); ++i) {
                    int init_val = initial[i];
                    if (init_val > 0) {
                        arcs.push_back(std::make_pair(Vector2i(-1, i), initial[i]));
                        supply += init_val;
                    } else if (init_val < 0) {
                        demand -= init_val;
                        arcs.push_back(std::make_pair(Vector2i(i, initial.size()), -init_val));
                    }
                }

                std::unique_ptr<MaxFlowHelper> solver = nullptr;
                if (supply < 20) {
                    solver = std::make_unique<ECMaxFlowHelper>();
                } else {
                    solver = std::make_unique<BoykovMaxFlowHelper>();
                }

                solver->resize(initial.size() + 2, arc_ids.size());

                std::set<int> ids;
                for (int i = 0; i < arcs.size(); ++i) {
                    int v1 = arcs[i].first[0] + 1;
                    int v2 = arcs[i].first[1] + 1;
                    int c = arcs[i].second;
                    if (v1 == 0 || v2 == initial.size() + 1) {
                        solver->addEdge(v1, v2, c, 0, -1);
                    } else {
                        if (arc_ids[i] > 0)
                            solver->addEdge(v1, v2, std::max(0, c + edge_capacity),
                                            std::max(0, -c + edge_capacity), arc_ids[i] - 1);
                        else {
                            if (c > 0)
                                solver->addEdge(v1, v2, std::max(0, c - 1),
                                                std::max(0, -c + edge_capacity), -1 - arc_ids[i]);
                            else
                                solver->addEdge(v1, v2, std::max(0, c + edge_capacity),
                                                std::max(0, -c - 1), -1 - arc_ids[i]);
                        }
                    }
                }
                int flow_count = solver->compute();

                solver->applyTo(EdgeDiff);

                //                lprintf("flow_count = %d, supply = %d\n", flow_count, supply);
                if (flow_count == supply) fullFlow = true;
                if (level != 0 || fullFlow) break;
                edge_capacity += 1;
                iter++;
                if (iter == 10) {
                    /* Probably won't converge. */
                    break;
                }
                //                lprintf("Not full flow, edge_capacity += 1\n");
            }

            if (level != 0) {
                auto &nEdgeDiff = mRes.mEdgeDiff[level - 1];
                auto &toUpper = mRes.mToUpperEdges[level - 1];
                auto &toUpperOrients = mRes.mToUpperOrients[level - 1];
                for (int i = 0; i < toUpper.size(); ++i) {
                    if (toUpper[i] >= 0) {
                        int orient = (4 - toUpperOrients[i]) % 4;
                        nEdgeDiff[i] = mathext::rshift90(EdgeDiff[toUpper[i]], orient);
                    }
                }
            }
        }
    }
}


namespace quadriflow {
    entities::Mesh to_mesh(const Parametrizer &field) {
        entities::Mesh mesh;

        for (int i = 0; i < field.m_vertices.cols(); ++i) {
            const auto vh = mesh.add_vertex(entities::Mesh::Point(
                field.m_vertices(0, i),
                field.m_vertices(1, i),
                field.m_vertices(2, i)
            ));

            assert(vh.idx() == i);
        }

        auto F = field.m_faces;
        int face_valence = F.rows();
        for (int f = 0; f < F.cols(); ++f) {
            const auto fh = mesh.add_face({
                entities::Mesh::VertexHandle(F(0, f)),
                entities::Mesh::VertexHandle(F(1, f)),
                entities::Mesh::VertexHandle(F(2, f))
            });

            assert(fh.idx() == f);

            //     for (auto i = 0; i < face_valence; ++i) {
            //         auto v = F(i, f);
            //         auto v_next = F((i + 1) % face_valence, f);
            //         auto e = face_valence * f + i;
            //
            //         const auto heh = mesh.halfedge_handle(mesh.vertex_handle(v));
            //         const auto v_to = mesh.to_vertex_handle(heh);
            //
            //         assert(heh.idx() == e);
            //         assert(v_to.idx() == v_next);
            //     }
        }

        return mesh;
    }

    void generate_adjacency_matrix_uniform(
        const MatrixXi &F,
        const VectorXi &V2E,
        const VectorXi &E2E,
        const VectorXi &nonManifold,
        AdjacentMatrix &adj
    ) {
        adj.resize(V2E.size());
        for (int i = 0; i < adj.size(); ++i) {
            int start = V2E[i];
            int edge = start;
            if (start == -1)
                continue;
            do {
                int base = edge % 3, f = edge / 3;
                int opp = E2E[edge], next = mathext::dedge_next_3(opp);
                if (adj[i].empty())
                    adj[i].push_back(Link(F((base + 2) % 3, f)));
                if (opp == -1 || next != start) {
                    adj[i].push_back(Link(F((base + 1) % 3, f)));
                    if (opp == -1)
                        break;
                }
                edge = next;
            } while (edge != start);
        }
    }

    void subdivide_edges_to_length(
        MatrixXi &F,
        MatrixXd &V,
        VectorXd &rho,
        VectorXi &V2E,
        VectorXi &E2E,
        VectorXi &boundary,
        VectorXi &nonmanifold,
        double maxLength
    ) {
        typedef std::pair<double, int> Edge;

        std::priority_queue<Edge> queue;

        maxLength *= maxLength;

        for (int i = 0; i < E2E.size(); ++i) {
            int v0 = F(i % 3, i / 3), v1 = F((i + 1) % 3, i / 3);
            if (nonmanifold[v0] || nonmanifold[v1]) continue;
            double length = (V.col(v0) - V.col(v1)).squaredNorm();
            if (length > maxLength || length > std::max(maxLength * 0.75, std::min(rho[v0], rho[v1]) * 1.0)) {
                int other = E2E[i];
                if (other == -1 || other > i) queue.push(Edge(length, i));
            }
        }

        spdlog::debug("Subdividing {} edges to length {}", queue.size(), maxLength);

        int nV = V.cols(), nF = F.cols(), nSplit = 0;
        /*
        /   v0  \
        v1p 1 | 0 v0p
        \   v1  /

        /   v0  \
        /  1 | 0  \
        v1p - vn - v0p
        \  2 | 3  /
        \   v1  /

        f0: vn, v0p, v0
        f1: vn, v0, v1p
        f2: vn, v1p, v1
        f3: vn, v1, v0p
        */
        int counter = 0;
        while (!queue.empty()) {
            counter += 1;
            Edge edge = queue.top();
            queue.pop();
            int e0 = edge.second, e1 = E2E[e0];
            bool is_boundary = e1 == -1;
            int f0 = e0 / 3, f1 = is_boundary ? -1 : (e1 / 3);
            int v0 = F(e0 % 3, f0), v0p = F((e0 + 2) % 3, f0), v1 = F((e0 + 1) % 3, f0);
            if ((V.col(v0) - V.col(v1)).squaredNorm() != edge.first) {
                continue;
            }
            int v1p = is_boundary ? -1 : F((e1 + 2) % 3, f1);
            int vn = nV++;
            nSplit++;
            /* Update m_vertices */
            if (nV > V.cols()) {
                V.conservativeResize(V.rows(), V.cols() * 2);
                rho.conservativeResize(V.cols() * 2);
                V2E.conservativeResize(V.cols());
                boundary.conservativeResize(V.cols());
                nonmanifold.conservativeResize(V.cols());
            }

            /* Update m_vertices */
            V.col(vn) = (V.col(v0) + V.col(v1)) * 0.5f;
            rho[vn] = 0.5f * (rho[v0], rho[v1]);
            nonmanifold[vn] = false;
            boundary[vn] = is_boundary;

            /* Update F and m_E2E */
            int f2 = is_boundary ? -1 : (nF++);
            int f3 = nF++;
            if (nF > F.cols()) {
                F.conservativeResize(F.rows(), std::max(nF, (int) F.cols() * 2));
                E2E.conservativeResize(F.cols() * 3);
            }

            /* Update F */
            F.col(f0) << vn, v0p, v0;
            if (!is_boundary) {
                F.col(f1) << vn, v0, v1p;
                F.col(f2) << vn, v1p, v1;
            }
            F.col(f3) << vn, v1, v0p;

            /* Update m_E2E */
            const int e0p = E2E[mathext::dedge_prev_3(e0)], e0n = E2E[mathext::dedge_next_3(e0)];

#define sE2E(a, b) \
    E2E[a] = b;    \
    if (b != -1) E2E[b] = a;
            sE2E(3 * f0 + 0, 3 * f3 + 2);
            sE2E(3 * f0 + 1, e0p);
            sE2E(3 * f3 + 1, e0n);
            if (is_boundary) {
                sE2E(3 * f0 + 2, -1);
                sE2E(3 * f3 + 0, -1);
            } else {
                const int e1p = E2E[mathext::dedge_prev_3(e1)], e1n = E2E[mathext::dedge_next_3(e1)];
                sE2E(3 * f0 + 2, 3 * f1 + 0);
                sE2E(3 * f1 + 1, e1n);
                sE2E(3 * f1 + 2, 3 * f2 + 0);
                sE2E(3 * f2 + 1, e1p);
                sE2E(3 * f2 + 2, 3 * f3 + 0);
            }
#undef sE2E

            /* Update m_V2E */
            V2E[v0] = 3 * f0 + 2;
            V2E[vn] = 3 * f0 + 0;
            V2E[v1] = 3 * f3 + 1;
            V2E[v0p] = 3 * f0 + 1;
            if (!is_boundary) V2E[v1p] = 3 * f1 + 2;

            auto schedule = [&](int f) {
                for (int i = 0; i < 3; ++i) {
                    double length = (V.col(F(i, f)) - V.col(F((i + 1) % 3, f))).squaredNorm();
                    if (length > maxLength
                        || length > std::max(maxLength * 0.75, std::min(rho[F(i, f)], rho[F((i + 1) % 3, f)]) * 1.0))
                        queue.push(Edge(length, f * 3 + i));
                }
            };

            schedule(f0);
            if (!is_boundary) {
                schedule(f2);
                schedule(f1);
            };
            schedule(f3);
        }
        F.conservativeResize(F.rows(), nF);
        V.conservativeResize(V.rows(), nV);
        rho.conservativeResize(nV);
        V2E.conservativeResize(nV);
        boundary.conservativeResize(nV);
        nonmanifold.conservativeResize(nV);
        E2E.conservativeResize(nF * 3);
    }


    void subdivide_edge_to_length_considering_edge_differences(
        MatrixXi &F,
        MatrixXd &V,
        MatrixXd &N,
        MatrixXd &Q,
        MatrixXd &O,
        MatrixXd *S,
        VectorXi &V2E,
        VectorXi &E2E,
        VectorXi &boundary,
        VectorXi &nonmanifold,
        std::vector<Vector2i> &edge_diff,
        std::vector<DEdge> &edge_values,
        std::vector<Vector3i> &face_edgeOrients,
        std::vector<Vector3i> &face_edgeIds,
        std::vector<int> &sharp_edges,
        std::map<int, int> &singularities,
        int max_len
    ) {
        struct EdgeLink {
            int id;
            double length;
            Vector2i diff;

            int maxlen() const { return std::max(abs(diff[0]), abs(diff[1])); }

            bool operator<(const EdgeLink &link) const { return maxlen() < link.maxlen(); }
        };

        struct FaceOrient {
            int orient;
            Vector3i d;
            Vector3d q;
            Vector3d n;
        };

        std::vector<FaceOrient> face_spaces(F.cols());
        std::priority_queue<EdgeLink> queue;
        std::vector<Vector2i> diffs(E2E.size());
        for (int i = 0; i < F.cols(); ++i) {
            for (int j = 0; j < 3; ++j) {
                int eid = i * 3 + j;
                diffs[eid] = mathext::rshift90(edge_diff[face_edgeIds[i][j]], face_edgeOrients[i][j]);
            }
        }
        for (int i = 0; i < F.cols(); ++i) {
            FaceOrient orient{};
            orient.q = Q.col(F(0, i));
            orient.n = N.col(F(0, i));
            int orient_diff[3];
            for (int j = 0; j < 3; ++j) {
                int final_orient = face_edgeOrients[i][j];
                int eid = face_edgeIds[i][j];
                auto value = mathext::compat_orientation_extrinsic_index_4(
                    Q.col(edge_values[eid].x), N.col(edge_values[eid].x), orient.q, orient.n);
                int target_orient = (value.second - value.first + 4) % 4;
                if (F(j, i) == edge_values[eid].y) target_orient = (target_orient + 2) % 4;
                orient_diff[j] = (final_orient - target_orient + 4) % 4;
            }
            if (orient_diff[0] == orient_diff[1])
                orient.orient = orient_diff[0];
            else if (orient_diff[0] == orient_diff[2])
                orient.orient = orient_diff[2];
            else if (orient_diff[1] == orient_diff[2])
                orient.orient = orient_diff[1];
            orient.d = Vector3i((orient_diff[0] - orient.orient + 4) % 4,
                                (orient_diff[1] - orient.orient + 4) % 4,
                                (orient_diff[2] - orient.orient + 4) % 4);
            face_spaces[i] = (orient);
        }
        for (int i = 0; i < E2E.size(); ++i) {
            int v0 = F(i % 3, i / 3), v1 = F((i + 1) % 3, i / 3);
            if (nonmanifold[v0] || nonmanifold[v1]) continue;
            double length = (V.col(v0) - V.col(v1)).squaredNorm();
            Vector2i diff = diffs[i];
            if (abs(diff[0]) > max_len || abs(diff[1]) > max_len) {
                int other = E2E[i];
                if (other == -1 || other > i) {
                    EdgeLink e;
                    e.id = i;
                    e.length = length;
                    e.diff = diff;
                    queue.push(e);
                }
            }
        }
        auto AnalyzeOrient = [&](int f0, const Vector3i &d) {
            for (int j = 0; j < 3; ++j) {
                int orient = face_spaces[f0].orient + d[j];
                int v = std::min(F(j, f0), F((j + 1) % 3, f0));
                auto value = mathext::compat_orientation_extrinsic_index_4(
                    Q.col(v), N.col(v), face_spaces[f0].q, face_spaces[f0].n);
                if (F(j, f0) != v) orient += 2;
                face_edgeOrients[f0][j] = (orient + value.second - value.first + 4) % 4;
            }
            face_spaces[f0].d = d;
            for (int j = 0; j < 3; ++j) {
                int eid = face_edgeIds[f0][j];
                int orient = face_edgeOrients[f0][j];
                auto diff = mathext::rshift90(diffs[f0 * 3 + j], (4 - orient) % 4);
                edge_diff[eid] = diff;
            }
        };
        auto FixOrient = [&](int f0) {
            for (int j = 0; j < 3; ++j) {
                auto diff = edge_diff[face_edgeIds[f0][j]];
                if (mathext::rshift90(diff, face_edgeOrients[f0][j]) != diffs[f0 * 3 + j]) {
                    int orient = 0;
                    while (orient < 4 && mathext::rshift90(diff, orient) != diffs[f0 * 3 + j]) orient += 1;
                    face_spaces[f0].d[j] =
                            (face_spaces[f0].d[j] + orient - face_edgeOrients[f0][j]) % 4;
                    face_edgeOrients[f0][j] = orient;
                }
            }
        };

        int nV = V.cols(), nF = F.cols(), nSplit = 0;
        /*
         /   v0  \
         v1p 1 | 0 v0p
         \   v1  /

         /   v0  \
         /  1 | 0  \
         v1p - vn - v0p
         \  2 | 3  /
         \   v1  /

         f0: vn, v0p, v0
         f1: vn, v0, v1p
         f2: vn, v1p, v1
         f3: vn, v1, v0p
         */
        int counter = 0;
        while (!queue.empty()) {
            counter += 1;
            EdgeLink edge = queue.top();
            queue.pop();

            int e0 = edge.id, e1 = E2E[e0];
            bool is_boundary = e1 == -1;
            int f0 = e0 / 3, f1 = is_boundary ? -1 : (e1 / 3);
            int v0 = F(e0 % 3, f0), v0p = F((e0 + 2) % 3, f0), v1 = F((e0 + 1) % 3, f0);
            if ((V.col(v0) - V.col(v1)).squaredNorm() != edge.length) {
                continue;
            }
            if (abs(diffs[e0][0]) < 2 && abs(diffs[e0][1]) < 2) continue;
            if (f1 != -1) {
                face_edgeOrients.push_back(Vector3i());
                sharp_edges.push_back(0);
                sharp_edges.push_back(0);
                sharp_edges.push_back(0);
                face_edgeIds.push_back(Vector3i());
            }
            int v1p = is_boundary ? -1 : F((e1 + 2) % 3, f1);
            int vn = nV++;
            nSplit++;
            if (nV > V.cols()) {
                V.conservativeResize(V.rows(), V.cols() * 2);
                N.conservativeResize(N.rows(), N.cols() * 2);
                Q.conservativeResize(Q.rows(), Q.cols() * 2);
                O.conservativeResize(O.rows(), O.cols() * 2);
                if (S)
                    S->conservativeResize(S->rows(), S->cols() * 2);
                V2E.conservativeResize(V.cols());
                boundary.conservativeResize(V.cols());
                nonmanifold.conservativeResize(V.cols());
            }

            V.col(vn) = (V.col(v0) + V.col(v1)) * 0.5;
            N.col(vn) = N.col(v0);
            Q.col(vn) = Q.col(v0);
            O.col(vn) = (O.col(v0) + O.col(v1)) * 0.5;
            if (S)
                S->col(vn) = S->col(v0);

            nonmanifold[vn] = false;
            boundary[vn] = is_boundary;

            int eid = face_edgeIds[f0][e0 % 3];
            int sharp_eid = sharp_edges[e0];
            int eid01 = face_edgeIds[f0][(e0 + 1) % 3];
            int sharp_eid01 = sharp_edges[f0 * 3 + (e0 + 1) % 3];
            int eid02 = face_edgeIds[f0][(e0 + 2) % 3];
            int sharp_eid02 = sharp_edges[f0 * 3 + (e0 + 2) % 3];

            int eid0, eid1, eid0p, eid1p;
            int sharp_eid0, sharp_eid1, sharp_eid0p, sharp_eid1p;

            eid0 = eid;
            sharp_eid0 = sharp_eid;
            edge_values[eid0] = DEdge(v0, vn);

            eid1 = edge_values.size();
            sharp_eid1 = sharp_eid;
            edge_values.push_back(DEdge(vn, v1));
            edge_diff.push_back(Vector2i());

            eid0p = edge_values.size();
            sharp_eid0p = 0;
            edge_values.push_back(DEdge(vn, v0p));
            edge_diff.push_back(Vector2i());

            int f2 = is_boundary ? -1 : (nF++);
            int f3 = nF++;
            sharp_edges.push_back(0);
            sharp_edges.push_back(0);
            sharp_edges.push_back(0);
            face_edgeIds.push_back(Vector3i());
            face_edgeOrients.push_back(Vector3i());

            if (nF > F.cols()) {
                F.conservativeResize(F.rows(), std::max(nF, (int) F.cols() * 2));
                face_spaces.resize(F.cols());
                E2E.conservativeResize(F.cols() * 3);
                diffs.resize(F.cols() * 3);
            }

            auto D01 = diffs[e0];
            auto D1p = diffs[e0 / 3 * 3 + (e0 + 1) % 3];
            auto Dp0 = diffs[e0 / 3 * 3 + (e0 + 2) % 3];

            Vector2i D0n = D01 / 2;

            auto orients1 = face_spaces[f0];
            F.col(f0) << vn, v0p, v0;
            face_edgeIds[f0] = Vector3i(eid0p, eid02, eid0);
            sharp_edges[f0 * 3] = sharp_eid0p;
            sharp_edges[f0 * 3 + 1] = sharp_eid02;
            sharp_edges[f0 * 3 + 2] = sharp_eid0;

            diffs[f0 * 3] = D01 + D1p - D0n;
            diffs[f0 * 3 + 1] = Dp0;
            diffs[f0 * 3 + 2] = D0n;
            int o1 = e0 % 3, o2 = e1 % 3;
            AnalyzeOrient(f0, Vector3i(0, orients1.d[(o1 + 2) % 3], orients1.d[o1]));
            if (!is_boundary) {
                auto orients2 = face_spaces[f1];
                int eid11 = face_edgeIds[f1][(e1 + 1) % 3];
                int sharp_eid11 = sharp_edges[f1 * 3 + (e1 + 1) % 3];
                int eid12 = face_edgeIds[f1][(e1 + 2) % 3];
                int sharp_eid12 = sharp_edges[f1 * 3 + (e1 + 2) % 3];

                auto Ds10 = diffs[e1];
                auto Ds0p = diffs[e1 / 3 * 3 + (e1 + 1) % 3];

                auto Dsp1 = diffs[e1 / 3 * 3 + (e1 + 2) % 3];
                int orient = 0;
                while (mathext::rshift90(D01, orient) != Ds10) orient += 1;
                Vector2i Dsn0 = mathext::rshift90(D0n, orient);

                F.col(f1) << vn, v0, v1p;
                eid1p = edge_values.size();
                sharp_eid1p = 0;
                edge_values.push_back(DEdge(vn, v1p));
                edge_diff.push_back(Vector2i());

                sharp_edges[f1 * 3] = sharp_eid0;
                sharp_edges[f1 * 3 + 1] = sharp_eid11;
                sharp_edges[f1 * 3 + 2] = sharp_eid1p;
                face_edgeIds[f1] = (Vector3i(eid0, eid11, eid1p));
                diffs[f1 * 3] = Dsn0;
                diffs[f1 * 3 + 1] = Ds0p;
                diffs[f1 * 3 + 2] = Dsp1 + (Ds10 - Dsn0);

                AnalyzeOrient(f1, Vector3i(orients2.d[o2], orients2.d[(o2 + 1) % 3], 0));

                face_spaces[f2] = face_spaces[f1];
                sharp_edges[f2 * 3] = sharp_eid1p;
                sharp_edges[f2 * 3 + 1] = sharp_eid12;
                sharp_edges[f2 * 3 + 2] = sharp_eid1;
                face_edgeIds[f2] = (Vector3i(eid1p, eid12, eid1));
                F.col(f2) << vn, v1p, v1;
                diffs[f2 * 3] = -Dsp1 - (Ds10 - Dsn0);
                diffs[f2 * 3 + 1] = Dsp1;
                diffs[f2 * 3 + 2] = Ds10 - Dsn0;

                AnalyzeOrient(f2, Vector3i(0, orients2.d[(o2 + 2) % 3], orients2.d[o2]));
            }
            face_spaces[f3] = face_spaces[f0];
            sharp_edges[f3 * 3] = sharp_eid1;
            sharp_edges[f3 * 3 + 1] = sharp_eid01;
            sharp_edges[f3 * 3 + 2] = sharp_eid0p;
            face_edgeIds[f3] = (Vector3i(eid1, eid01, eid0p));
            F.col(f3) << vn, v1, v0p;
            diffs[f3 * 3] = D01 - D0n;
            diffs[f3 * 3 + 1] = D1p;
            diffs[f3 * 3 + 2] = D0n - (D01 + D1p);

            AnalyzeOrient(f3, Vector3i(orients1.d[o1], orients1.d[(o1 + 1) % 3], 0));

            FixOrient(f0);
            if (!is_boundary) {
                FixOrient(f1);
                FixOrient(f2);
            }
            FixOrient(f3);

            const int e0p = E2E[mathext::dedge_prev_3(e0)], e0n = E2E[mathext::dedge_next_3(e0)];

#define sE2E(a, b) \
    E2E[a] = b;    \
    if (b != -1) E2E[b] = a;
            sE2E(3 * f0 + 0, 3 * f3 + 2);
            sE2E(3 * f0 + 1, e0p);
            sE2E(3 * f3 + 1, e0n);
            if (is_boundary) {
                sE2E(3 * f0 + 2, -1);
                sE2E(3 * f3 + 0, -1);
            } else {
                const int e1p = E2E[mathext::dedge_prev_3(e1)], e1n = E2E[mathext::dedge_next_3(e1)];
                sE2E(3 * f0 + 2, 3 * f1 + 0);
                sE2E(3 * f1 + 1, e1n);
                sE2E(3 * f1 + 2, 3 * f2 + 0);
                sE2E(3 * f2 + 1, e1p);
                sE2E(3 * f2 + 2, 3 * f3 + 0);
            }
#undef sE2E

            V2E[v0] = 3 * f0 + 2;
            V2E[vn] = 3 * f0 + 0;
            V2E[v1] = 3 * f3 + 1;
            V2E[v0p] = 3 * f0 + 1;
            if (!is_boundary) V2E[v1p] = 3 * f1 + 2;

            auto schedule = [&](int f) {
                for (int i = 0; i < 3; ++i) {
                    if (abs(diffs[f * 3 + i][0]) > max_len || abs(diffs[f * 3 + i][1]) > max_len) {
                        EdgeLink e;
                        e.id = f * 3 + i;
                        e.length = (V.col(F((i + 1) % 3, f)) - V.col(F(i, f))).squaredNorm();
                        e.diff = diffs[f * 3 + i];
                        queue.push(e);
                    }
                }
            };

            schedule(f0);
            if (!is_boundary) {
                schedule(f2);
                schedule(f1);
            };
            schedule(f3);
        }
        F.conservativeResize(F.rows(), nF);
        V.conservativeResize(V.rows(), nV);
        N.conservativeResize(V.rows(), nV);
        Q.conservativeResize(V.rows(), nV);
        O.conservativeResize(V.rows(), nV);
        if (S) {
            S->conservativeResize(S->rows(), nV);
        }
        V2E.conservativeResize(nV);
        boundary.conservativeResize(nV);
        nonmanifold.conservativeResize(nV);
        E2E.conservativeResize(nF * 3);
        for (int i = 0; i < F.cols(); ++i) {
            for (int j = 0; j < 3; ++j) {
                auto diff = edge_diff[face_edgeIds[i][j]];
                if (abs(diff[0]) > 1 || abs(diff[1]) > 1) {
                    printf("wrong init %d %d!\n", face_edgeIds[i][j], i * 3 + j);
                    exit(0);
                }
            }
        }
        for (int i = 0; i < edge_diff.size(); ++i) {
            if (abs(edge_diff[i][0]) > 1 || abs(edge_diff[i][1]) > 1) {
                printf("wrong...\n");
                exit(0);
            }
        }
    }

    void Parametrizer::initialize(
        entities::SDFn sdfn,
        int target_face_count,
        bool with_scale
    ) {
        m_hierarchy.clearConstraints();

        const auto [vertices_normalized, scale, offset] = mathext::normalize(m_vertices.cast<float>().transpose());
        m_vertices = vertices_normalized.cast<double>().transpose();
        m_normalize_scale = scale;
        m_normalize_offset = offset.cast<double>();

#ifdef DEV_DEBUG
        const entities::SDFn sdfn_scale = sdfn::scale(sdfn, scale, offset);
        const auto sdf = sdfn_scale(m_vertices.transpose().cast<float>());
        entities::Mesh mesh;
        for (int i = 0; i < m_vertices.cols(); ++i) {
            if (sdf[i] < 0) {
                mesh.add_vertex(entities::Mesh::Point(
                    m_vertices(0, i),
                    m_vertices(1, i),
                    m_vertices(2, i)
                ));
            }
        }

        OpenMesh::IO::write_mesh(mesh, "../tests/out/stage/6-normalized.ply");
#endif

        analyze_mesh();

        // initialize m_rho
        m_rho.resize(m_vertices.cols(), 1);
        m_rho.setConstant(1);

        // initialize the scale of the mesh
        if (target_face_count <= 0) {
            m_scale = sqrt(m_surface_area / m_vertices.cols());
        } else {
            m_scale = std::sqrt(m_surface_area / target_face_count);
        }

        // Computes the directed graph and subdivides if the scale is larger than the maximum edge length.
        double target_len = std::min(m_scale / 2, m_average_edge_length * 2);
        if (target_len < m_max_edge_length) {
            while (!mathext::compute_direct_graph(m_vertices, m_faces, m_V2E, m_E2E, m_boundary, m_non_manifold));
            subdivide_edges_to_length(m_faces, m_vertices, m_rho, m_V2E, m_E2E, m_boundary, m_non_manifold, target_len);
        }
        while (!mathext::compute_direct_graph(m_vertices, m_faces, m_V2E, m_E2E, m_boundary, m_non_manifold));

        // Compute the adjacency matrix
        generate_adjacency_matrix_uniform(m_faces, m_V2E, m_E2E, m_non_manifold, m_adjacency_matrix);

        // Computes the shortest edge per vertex. FIXME
        for (int iter = 0; iter < 5; ++iter) {
            VectorXd r(m_rho.size());
            for (int i = 0; i < m_rho.size(); ++i) {
                r[i] = m_rho[i];
                for (auto &id: m_adjacency_matrix[i]) {
                    r[i] = std::min(r[i], m_rho[id.id]);
                }
            }
            m_rho = r;
        }

        m_edges_preserve.resize(m_faces.cols() * 3);
        find_and_set_edges();
        compute_normals();
        compute_vertex_area();

        if (with_scale) {
            m_triangle_space.resize(m_faces.cols());
            for (int i = 0; i < m_faces.cols(); ++i) {
                Matrix3d p, q;
                p.col(0) = m_vertices.col(m_faces(1, i)) - m_vertices.col(m_faces(0, i));
                p.col(1) = m_vertices.col(m_faces(2, i)) - m_vertices.col(m_faces(0, i));
                p.col(2) = m_faces_normals.col(i);
                q = p.inverse();
                m_triangle_space[i].resize(2, 3);
                for (int j = 0; j < 2; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        m_triangle_space[i](j, k) = q(j, k);
                    }
                }
            }
        }

        m_hierarchy.m_vertex_area[0] = std::move(m_vertex_area);
        m_hierarchy.m_adjacency[0] = std::move(m_adjacency_matrix);
        m_hierarchy.m_normals[0] = std::move(m_normals_vertices);
        m_hierarchy.m_vertices[0] = std::move(m_vertices);
        m_hierarchy.m_E2E = std::move(m_E2E);
        m_hierarchy.m_faces = std::move(m_faces);
        m_hierarchy.Initialize(m_scale, with_scale);
    }

    void Parametrizer::analyze_mesh() {
        m_surface_area = 0;
        m_average_edge_length = 0;
        m_max_edge_length = 0;
        for (int f = 0; f < m_faces.cols(); ++f) {
            Vector3d v[3] = {
                m_vertices.col(m_faces(0, f)), m_vertices.col(m_faces(1, f)),
                m_vertices.col(m_faces(2, f))
            };
            double area = 0.5f * (v[1] - v[0]).cross(v[2] - v[0]).norm();
            m_surface_area += area;
            for (int i = 0; i < 3; ++i) {
                double len = (v[(i + 1) % 3] - v[i]).norm();
                m_average_edge_length += len;
                if (len > m_max_edge_length) m_max_edge_length = len;
            }
        }
        m_average_edge_length /= (m_faces.cols() * 3);
    }

    void Parametrizer::compute_vertex_area() {
        m_vertex_area.resize(m_vertices.cols());
        m_vertex_area.setZero();

        for (int i = 0; i < m_V2E.size(); ++i) {
            int edge = m_V2E[i], stop = edge;
            if (m_non_manifold[i] || edge == -1) continue;
            double vertex_area = 0;
            do {
                int ep = mathext::dedge_prev_3(edge), en = mathext::dedge_next_3(edge);

                Vector3d v = m_vertices.col(m_faces(edge % 3, edge / 3));
                Vector3d vn = m_vertices.col(m_faces(en % 3, en / 3));
                Vector3d vp = m_vertices.col(m_faces(ep % 3, ep / 3));

                Vector3d face_center = (v + vp + vn) * (1.0f / 3.0f);
                Vector3d prev = (v + vp) * 0.5f;
                Vector3d next = (v + vn) * 0.5f;

                vertex_area += 0.5f * ((v - prev).cross(v - face_center).norm() +
                                       (v - next).cross(v - face_center).norm());

                int opp = m_E2E[edge];
                if (opp == -1) break;
                edge = mathext::dedge_next_3(opp);
            } while (edge != stop);

            m_vertex_area[i] = vertex_area;
        }
    }

    void Parametrizer::compute_normals() {
        /* Compute face normals */
        m_faces_normals.resize(3, m_faces.cols());
        for (int f = 0; f < m_faces.cols(); ++f) {
            Vector3d v0 = m_vertices.col(m_faces(0, f)), v1 = m_vertices.col(m_faces(1, f)), v2 = m_vertices.col(
                        m_faces(2, f)),
                    n = (v1 - v0).cross(v2 - v0);
            double norm = n.norm();
            if (norm < RCPOVERFLOW) {
                n = Vector3d::UnitX();
            } else {
                n /= norm;
            }
            m_faces_normals.col(f) = n;
        }

        m_normals_vertices.resize(3, m_vertices.cols());
        for (int i = 0; i < m_V2E.rows(); ++i) {
            int edge = m_V2E[i];
            if (m_non_manifold[i] || edge == -1) {
                m_normals_vertices.col(i) = Vector3d::UnitX();
                continue;
            }

            int stop = edge;
            do {
                if (m_edges_preserve[edge]) break;
                edge = m_E2E[edge];
                if (edge != -1) edge = mathext::dedge_next_3(edge);
            } while (edge != stop && edge != -1);
            if (edge == -1)
                edge = stop;
            else
                stop = edge;
            Vector3d normal = Vector3d::Zero();
            do {
                int idx = edge % 3;

                Vector3d d0 = m_vertices.col(m_faces((idx + 1) % 3, edge / 3)) - m_vertices.col(i);
                Vector3d d1 = m_vertices.col(m_faces((idx + 2) % 3, edge / 3)) - m_vertices.col(i);
                double angle = mathext::fast_acos(d0.dot(d1) / std::sqrt(d0.squaredNorm() * d1.squaredNorm()));

                /* "Computing Vertex Normals from Polygonal Facets"
                 by Grit Thuermer and Charles A. Wuethrich, JGT 1998, Vol 3 */
                if (std::isfinite(angle)) normal += m_faces_normals.col(edge / 3) * angle;

                int opp = m_E2E[edge];
                if (opp == -1) break;

                edge = mathext::dedge_next_3(opp);
                if (m_edges_preserve[edge]) break;
            } while (edge != stop);
            double norm = normal.norm();
            m_normals_vertices.col(i) = norm > RCPOVERFLOW ? Vector3d(normal / norm) : Vector3d::UnitX();
        }
    }

    void Parametrizer::find_and_set_edges(
        entities::SDFn sdfn,
        const float max_angle
    ) {
        spdlog::debug("Finding edges via SDFn to preserve them in the optimization");

        entities::Mesh mesh = to_mesh(*this);
        const auto field_angular = smoothing::laplacian_angular_field(sdfn, mesh);

        int count = 0;
        for (int idx_v = 0; idx_v < field_angular.size(); ++idx_v) {
            if (field_angular[idx_v] >= max_angle) {
                for (auto vv = mesh.vv_begin(mesh.vertex_handle(idx_v));
                     vv != mesh.vv_end(mesh.vertex_handle(idx_v)); ++vv) {
                    if (field_angular[vv->idx()] >= max_angle) {
                        const auto edge_start = m_V2E[idx_v];
                        const auto edge_inverse = m_V2E[vv->idx()];
                        if (edge_start == m_E2E[edge_inverse]) {
                            count++;
                        }
                        // if (m_edges_preserve[e] == 0 && e == eb) {
                        // m_edges_preserve[e] = 1;
                        // count++;
                        // }
                    }
                }
            }
        }

#ifdef DEV_DEBUG
        std::vector<entities::Mesh::VertexHandle> to_delete;

        entities::Mesh mesh_edges;
        for (auto vh = mesh.vertices_begin(); vh != mesh.vertices_end(); ++vh) {
            if (field_angular[vh->idx()] >= max_angle) {
                mesh_edges.add_vertex(mesh.point(*vh));
            }
        }

        OpenMesh::IO::write_mesh(mesh_edges, "../tests/out/8-edges.ply");
        spdlog::debug("Found {} edges to preserve", count);
#endif
    }

    void Parametrizer::find_and_set_edges() {
        spdlog::debug("Finding boundaries to preserve them in the optimization");

        std::vector<Vector3d> face_normals(m_faces.cols());
        for (int i = 0; i < m_faces.cols(); ++i) {
            Vector3d p1 = m_vertices.col(m_faces(0, i));
            Vector3d p2 = m_vertices.col(m_faces(1, i));
            Vector3d p3 = m_vertices.col(m_faces(2, i));
            face_normals[i] = (p2 - p1).cross(p3 - p1).normalized();
        }

        const double threshold = cos(60.0 / 180.0 * 3.141592654);

        int count = 0;
        for (int i = 0; i < m_edges_preserve.size(); ++i) {
            Vector3d &n1 = face_normals[i / 3];
            Vector3d &n2 = face_normals[m_E2E[i] / 3];
            if (n1.dot(n2) < threshold) {
                m_edges_preserve[i] = 1;
                count++;
            }
        }

#ifdef DEV_DEBUG
        spdlog::debug("Found {} edges to preserve", count);
#endif
    }

    void Parametrizer::build_edge_info() {
        auto &F = m_hierarchy.m_faces;
        auto &E2E = m_hierarchy.m_E2E;

        m_edge_difference.clear();
        m_edge_values.clear();
        m_face_edge_ids.resize(F.cols(), Vector3i(-1, -1, -1));
        for (int i = 0; i < F.cols(); ++i) {
            for (int j = 0; j < 3; ++j) {
                int k1 = j, k2 = (j + 1) % 3;
                int v1 = F(k1, i);
                int v2 = F(k2, i);
                DEdge e2(v1, v2);
                Vector2i diff2;
                int rank2;
                if (v1 > v2) {
                    rank2 = m_singularity_rank(k2, i);
                    diff2 = mathext::rshift90(
                        Vector2i(
                            -m_singularity_index(k1 * 2, i),
                            -m_singularity_index(k1 * 2 + 1, i)
                        ),
                        rank2
                    );
                } else {
                    rank2 = m_singularity_rank(k1, i);
                    diff2 = mathext::rshift90(
                        Vector2i(
                            m_singularity_index(k1 * 2, i),
                            m_singularity_index(k1 * 2 + 1, i)
                        ),
                        rank2
                    );
                }
                int current_eid = i * 3 + k1;
                int eid = E2E[current_eid];
                int eID1 = m_face_edge_ids[current_eid / 3][current_eid % 3];
                int eID2 = -1;
                if (eID1 == -1) {
                    eID2 = m_edge_values.size();
                    m_edge_values.push_back(e2);
                    m_edge_difference.push_back(diff2);
                    m_face_edge_ids[i][k1] = eID2;
                    if (eid != -1) m_face_edge_ids[eid / 3][eid % 3] = eID2;
                } else if (!m_singularities.count(i)) {
                    eID2 = m_face_edge_ids[eid / 3][eid % 3];
                    m_edge_difference[eID2] = diff2;
                }
            }
        }
    }

    void Parametrizer::build_integer_constraints() {
        auto &F = m_hierarchy.m_faces;
        auto &Q = m_hierarchy.m_orientation[0];
        auto &N = m_hierarchy.m_normals[0];
        m_face_edge_orientation.resize(F.cols());

        //Random number generator (for shuffling)
        std::random_device rd;
        std::mt19937 g(rd());
        g.seed(m_hierarchy.rng_seed);

        // undirected edge to direct edge
        std::vector<std::pair<int, int> > E2D(m_edge_difference.size(), std::make_pair(-1, -1));
        for (int i = 0; i < F.cols(); ++i) {
            int v0 = F(0, i);
            int v1 = F(1, i);
            int v2 = F(2, i);
            DEdge e0(v0, v1), e1(v1, v2), e2(v2, v0);
            const Vector3i &eid = m_face_edge_ids[i];
            Vector2i variable_id[3];
            for (int i = 0; i < 3; ++i) {
                variable_id[i] = Vector2i(eid[i] * 2 + 1, eid[i] * 2 + 2);
            }
            auto index1 =
                    mathext::compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v1), N.col(v1));
            auto index2 =
                    mathext::compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v2), N.col(v2));

            int rank1 = (index1.first - index1.second + 4) % 4; // v1 -> v0
            int rank2 = (index2.first - index2.second + 4) % 4; // v2 -> v0
            int orients[3] = {0}; // == {0, 0, 0}
            if (v1 < v0) {
                variable_id[0] = -mathext::rshift90(variable_id[0], rank1);
                orients[0] = (rank1 + 2) % 4;
            } else {
                orients[0] = 0;
            }
            if (v2 < v1) {
                variable_id[1] = -mathext::rshift90(variable_id[1], rank2);
                orients[1] = (rank2 + 2) % 4;
            } else {
                variable_id[1] = mathext::rshift90(variable_id[1], rank1);
                orients[1] = rank1;
            }
            if (v2 < v0) {
                variable_id[2] = mathext::rshift90(variable_id[2], rank2);
                orients[2] = rank2;
            } else {
                variable_id[2] = -variable_id[2];
                orients[2] = 2;
            }
            m_face_edge_orientation[i] = Vector3i(orients[0], orients[1], orients[2]);
            for (int j = 0; j < 3; ++j) {
                int eid = m_face_edge_ids[i][j];
                if (E2D[eid].first == -1)
                    E2D[eid].first = i * 3 + j;
                else
                    E2D[eid].second = i * 3 + j;
            }
        }

        // a face disajoint tree
        DisjointOrientTree disajoint_orient_tree = DisjointOrientTree(F.cols());
        // merge the whole face graph except for the singularity in which there exists a spanning tree
        // which contains consistent orientation
        std::vector<int> sharpUE(E2D.size());
        for (int i = 0; i < m_edges_preserve.size(); ++i) {
            if (m_edges_preserve[i]) {
                sharpUE[m_face_edge_ids[i / 3][i % 3]] = 1;
            }
        }

        for (int i = 0; i < E2D.size(); ++i) {
            auto &edge_c = E2D[i];
            int f0 = edge_c.first / 3;
            int f1 = edge_c.second / 3;
            if (edge_c.first == -1 || edge_c.second == -1) continue;
            if (m_singularities.count(f0) || m_singularities.count(f1) || sharpUE[i]) continue;
            int orient1 = m_face_edge_orientation[f0][edge_c.first % 3];
            int orient0 = (m_face_edge_orientation[f1][edge_c.second % 3] + 2) % 4;
            disajoint_orient_tree.Merge(f0, f1, orient0, orient1);
        }

        // merge singularity later
        for (auto &f: m_singularities) {
            for (int i = 0; i < 3; ++i) {
                if (sharpUE[m_face_edge_ids[f.first][i]]) continue;
                auto &edge_c = E2D[m_face_edge_ids[f.first][i]];
                if (edge_c.first == -1 || edge_c.second == -1) continue;
                int v0 = edge_c.first / 3;
                int v1 = edge_c.second / 3;
                int orient1 = m_face_edge_orientation[v0][edge_c.first % 3];
                int orient0 = (m_face_edge_orientation[v1][edge_c.second % 3] + 2) % 4;
                disajoint_orient_tree.Merge(v0, v1, orient0, orient1);
            }
        }

        for (int i = 0; i < sharpUE.size(); ++i) {
            if (sharpUE[i] == 0) continue;
            auto &edge_c = E2D[i];
            if (edge_c.first == -1 || edge_c.second == -1) continue;
            int f0 = edge_c.first / 3;
            int f1 = edge_c.second / 3;
            int orient1 = m_face_edge_orientation[f0][edge_c.first % 3];
            int orient0 = (m_face_edge_orientation[f1][edge_c.second % 3] + 2) % 4;
            disajoint_orient_tree.Merge(f0, f1, orient0, orient1);
        }

        // all the face has the same parent.  we rotate every face to the space of that parent.
        for (int i = 0; i < m_face_edge_orientation.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                m_face_edge_orientation[i][j] =
                        (m_face_edge_orientation[i][j] + disajoint_orient_tree.Orient(i)) % 4;
            }
        }

        std::vector<int> sharp_colors(m_face_edge_ids.size(), -1);
        int num_sharp_component = 0;
        // label the connected component connected by non-fixed edges
        // we need this because we need sink flow (demand) == source flow (supply) for each component
        // rather than global
        for (int i = 0; i < sharp_colors.size(); ++i) {
            if (sharp_colors[i] != -1) continue;
            sharp_colors[i] = num_sharp_component;
            std::queue<int> q;
            q.push(i);
            int counter = 0;
            while (!q.empty()) {
                int v = q.front();
                q.pop();
                for (int i = 0; i < 3; ++i) {
                    int e = m_face_edge_ids[v][i];
                    int deid1 = E2D[e].first;
                    int deid2 = E2D[e].second;
                    if (deid1 == -1 || deid2 == -1) continue;
                    if (abs(m_face_edge_orientation[deid1 / 3][deid1 % 3] -
                            m_face_edge_orientation[deid2 / 3][deid2 % 3] + 4) %
                        4 !=
                        2 ||
                        sharpUE[e]) {
                        continue;
                    }
                    for (int k = 0; k < 2; ++k) {
                        int f = (k == 0) ? E2D[e].first / 3 : E2D[e].second / 3;
                        if (sharp_colors[f] == -1) {
                            sharp_colors[f] = num_sharp_component;
                            q.push(f);
                        }
                    }
                }
                counter += 1;
            }
            num_sharp_component += 1;
        } {
            std::vector<int> total_flows(num_sharp_component);
            // check if each component is full-flow
            for (int i = 0; i < m_face_edge_ids.size(); ++i) {
                Vector2i diff(0, 0);
                for (int j = 0; j < 3; ++j) {
                    int orient = m_face_edge_orientation[i][j];
                    diff += mathext::rshift90(m_edge_difference[m_face_edge_ids[i][j]], orient);
                }
                total_flows[sharp_colors[i]] += diff[0] + diff[1];
            }

            // build "variable"
            m_variables.resize(m_edge_difference.size() * 2, std::make_pair(Vector2i(-1, -1), 0));
            for (int i = 0; i < m_face_edge_ids.size(); ++i) {
                for (int j = 0; j < 3; ++j) {
                    Vector2i sign = mathext::rshift90(Vector2i(1, 1), m_face_edge_orientation[i][j]);
                    int eid = m_face_edge_ids[i][j];
                    Vector2i index = mathext::rshift90(Vector2i(eid * 2, eid * 2 + 1), m_face_edge_orientation[i][j]);
                    for (int k = 0; k < 2; ++k) {
                        auto &p = m_variables[abs(index[k])];
                        if (p.first[0] == -1)
                            p.first[0] = i * 2 + k;
                        else
                            p.first[1] = i * 2 + k;
                        p.second += sign[k];
                    }
                }
            }

            // fixed variable that might be manually modified.
            // modified_variables[component_od][].first = fixed_variable_id
            // modified_variables[component_od][].second = 1 if two positive signs -1 if two negative
            // signs
            std::vector<std::vector<std::pair<int, int> > > modified_variables[2];
            for (int i = 0; i < 2; ++i) modified_variables[i].resize(total_flows.size());
            for (int i = 0; i < m_variables.size(); ++i) {
                if ((m_variables[i].first[1] == -1 || m_variables[i].second != 0) &&
                    m_allow_changes[i] == 1) {
                    int find = sharp_colors[m_variables[i].first[0] / 2];
                    int step = std::abs(m_variables[i].second) % 2;
                    if (total_flows[find] > 0) {
                        if (m_variables[i].second > 0 && m_edge_difference[i / 2][i % 2] > -1) {
                            modified_variables[step][find].push_back(std::make_pair(i, -1));
                        }
                        if (m_variables[i].second < 0 && m_edge_difference[i / 2][i % 2] < 1) {
                            modified_variables[step][find].push_back(std::make_pair(i, 1));
                        }
                    } else if (total_flows[find] < 0) {
                        if (m_variables[i].second < 0 && m_edge_difference[i / 2][i % 2] > -1) {
                            modified_variables[step][find].push_back(std::make_pair(i, -1));
                        }
                        if (m_variables[i].second > 0 && m_edge_difference[i / 2][i % 2] < 1) {
                            modified_variables[step][find].push_back(std::make_pair(i, 1));
                        }
                    }
                }
            }

            // uniformly random manually modify variables so that the network has full flow.
            for (int i = 0; i < 2; ++i) {
                for (auto &modified_var: modified_variables[i]) {
                    std::shuffle(modified_var.begin(), modified_var.end(), g);
                }
            }

            for (int j = 0; j < total_flows.size(); ++j) {
                for (int ii = 0; ii < 2; ++ii) {
                    if (total_flows[j] == 0) continue;
                    int max_num;
                    if (ii == 0)
                        max_num =
                                std::min(abs(total_flows[j]) / 2, (int) modified_variables[ii][j].size());
                    else
                        max_num = std::min(abs(total_flows[j]), (int) modified_variables[ii][j].size());
                    int dir = (total_flows[j] > 0) ? -1 : 1;
                    for (int i = 0; i < max_num; ++i) {
                        auto &info = modified_variables[ii][j][i];
                        m_edge_difference[info.first / 2][info.first % 2] += info.second;
                        if (ii == 0)
                            total_flows[j] += 2 * dir;
                        else
                            total_flows[j] += dir;
                    }
                }
            }
        }

        std::vector<Vector4i> edge_to_constraints(E2D.size() * 2, Vector4i(-1, 0, -1, 0));
        for (int i = 0; i < m_face_edge_ids.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                int e = m_face_edge_ids[i][j];
                Vector2i index = mathext::rshift90(Vector2i(e * 2 + 1, e * 2 + 2), m_face_edge_orientation[i][j]);
                for (int k = 0; k < 2; ++k) {
                    int l = abs(index[k]);
                    int s = index[k] / l;
                    int ind = l - 1;
                    int equationID = i * 2 + k;
                    if (edge_to_constraints[ind][0] == -1) {
                        edge_to_constraints[ind][0] = equationID;
                        edge_to_constraints[ind][1] = s;
                    } else {
                        edge_to_constraints[ind][2] = equationID;
                        edge_to_constraints[ind][3] = s;
                    }
                }
            }
        }
        std::vector<std::pair<Vector2i, int> > arcs;
        std::vector<int> arc_ids;
        DisjointTree tree(m_face_edge_ids.size() * 2);
        for (int i = 0; i < edge_to_constraints.size(); ++i) {
            if (m_allow_changes[i] == 0) continue;
            if (edge_to_constraints[i][0] == -1 || edge_to_constraints[i][2] == -1) continue;
            if (edge_to_constraints[i][1] == -edge_to_constraints[i][3]) {
                int v1 = edge_to_constraints[i][0];
                int v2 = edge_to_constraints[i][2];
                tree.Merge(v1, v2);
                if (edge_to_constraints[i][1] < 0) std::swap(v1, v2);
                int current_v = m_edge_difference[i / 2][i % 2];
                arcs.push_back(std::make_pair(Vector2i(v1, v2), current_v));
            }
        }
        tree.BuildCompactParent();
        std::vector<int> total_flows(tree.CompactNum());
        // check if each component is full-flow
        for (int i = 0; i < m_face_edge_ids.size(); ++i) {
            Vector2i diff(0, 0);
            for (int j = 0; j < 3; ++j) {
                int orient = m_face_edge_orientation[i][j];
                diff += mathext::rshift90(m_edge_difference[m_face_edge_ids[i][j]], orient);
            }
            for (int j = 0; j < 2; ++j) {
                total_flows[tree.Index(i * 2 + j)] += diff[j];
            }
        }

        // build "variable"
        m_variables.resize(m_edge_difference.size() * 2);
        for (int i = 0; i < m_variables.size(); ++i) {
            m_variables[i].first = Vector2i(-1, -1);
            m_variables[i].second = 0;
        }
        for (int i = 0; i < m_face_edge_ids.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                Vector2i sign = mathext::rshift90(Vector2i(1, 1), m_face_edge_orientation[i][j]);
                int eid = m_face_edge_ids[i][j];
                Vector2i index = mathext::rshift90(Vector2i(eid * 2, eid * 2 + 1), m_face_edge_orientation[i][j]);
                for (int k = 0; k < 2; ++k) {
                    auto &p = m_variables[abs(index[k])];
                    if (p.first[0] == -1)
                        p.first[0] = i * 2 + k;
                    else
                        p.first[1] = i * 2 + k;
                    p.second += sign[k];
                }
            }
        }

        // fixed variable that might be manually modified.
        // modified_variables[component_od][].first = fixed_variable_id
        // modified_variables[component_od][].second = 1 if two positive signs -1 if two negative signs
        std::vector<std::vector<std::pair<int, int> > > modified_variables[2];
        for (int i = 0; i < 2; ++i) {
            modified_variables[i].resize(total_flows.size());
        }
        for (int i = 0; i < m_variables.size(); ++i) {
            if ((m_variables[i].first[1] == -1 || m_variables[i].second != 0) && m_allow_changes[i] == 1) {
                int find = tree.Index(m_variables[i].first[0]);
                int step = abs(m_variables[i].second) % 2;
                if (total_flows[find] > 0) {
                    if (m_variables[i].second > 0 && m_edge_difference[i / 2][i % 2] > -1) {
                        modified_variables[step][find].push_back(std::make_pair(i, -1));
                    }
                    if (m_variables[i].second < 0 && m_edge_difference[i / 2][i % 2] < 1) {
                        modified_variables[step][find].push_back(std::make_pair(i, 1));
                    }
                } else if (total_flows[find] < 0) {
                    if (m_variables[i].second < 0 && m_edge_difference[i / 2][i % 2] > -1) {
                        modified_variables[step][find].push_back(std::make_pair(i, -1));
                    }
                    if (m_variables[i].second > 0 && m_edge_difference[i / 2][i % 2] < 1) {
                        modified_variables[step][find].push_back(std::make_pair(i, 1));
                    }
                }
            }
        }

        // uniformly random manually modify variables so that the network has full flow.
        for (int j = 0; j < 2; ++j) {
            for (auto &modified_var: modified_variables[j])
                std::shuffle(modified_var.begin(), modified_var.end(), g);
        }
        for (int j = 0; j < total_flows.size(); ++j) {
            for (int ii = 0; ii < 2; ++ii) {
                if (total_flows[j] == 0) continue;
                int max_num;
                if (ii == 0)
                    max_num = std::min(abs(total_flows[j]) / 2, (int) modified_variables[ii][j].size());
                else
                    max_num = std::min(abs(total_flows[j]), (int) modified_variables[ii][j].size());
                int dir = (total_flows[j] > 0) ? -1 : 1;
                for (int i = 0; i < max_num; ++i) {
                    auto &info = modified_variables[ii][j][i];
                    m_edge_difference[info.first / 2][info.first % 2] += info.second;
                    if (ii == 0)
                        total_flows[j] += 2 * dir;
                    else
                        total_flows[j] += dir;
                }
            }
        }
    }

    void Parametrizer::compute_index_map(Hierarchy &hierarchy, int with_scale) {
        auto &V = hierarchy.m_vertices[0];
        auto &F = hierarchy.m_faces;
        auto &Q = hierarchy.m_orientation[0];
        auto &N = hierarchy.m_normals[0];
        auto &O = hierarchy.m_positions[0];
        auto &S = hierarchy.m_scales[0];

        build_edge_info();

        // Constraints for the integer optimization

        for (int i = 0; i < m_edges_preserve.size(); ++i) {
            if (m_edges_preserve[i]) {
                int e = m_face_edge_ids[i / 3][i % 3];
                if (m_edge_difference[e][0] * m_edge_difference[e][1] != 0) {
                    Vector3d d = O.col(m_edge_values[e].y) - O.col(m_edge_values[e].x);
                    Vector3d q = Q.col(m_edge_values[e].x);
                    Vector3d n = N.col(m_edge_values[e].x);
                    Vector3d qy = n.cross(q);
                    if (abs(q.dot(d)) > qy.dot(d)) {
                        m_edge_difference[e][1] = 0;
                    } else {
                        m_edge_difference[e][0] = 0;
                    }
                }
            }
        }
        std::map<int, std::pair<Vector3d, Vector3d> > sharp_constraints;
        std::set<int> sharpvert;
        for (int i = 0; i < m_edges_preserve.size(); ++i) {
            if (m_edges_preserve[i]) {
                sharpvert.insert(F(i % 3, i / 3));
                sharpvert.insert(F((i + 1) % 3, i / 3));
            }
        }

        m_allow_changes.resize(m_edge_difference.size() * 2, 1);
        for (int i = 0; i < m_edges_preserve.size(); ++i) {
            int e = m_face_edge_ids[i / 3][i % 3];
            if (sharpvert.count(m_edge_values[e].x) && sharpvert.count(m_edge_values[e].y)) {
                if (m_edges_preserve[i] != 0) {
                    for (int k = 0; k < 2; ++k) {
                        if (m_edge_difference[e][k] == 0) {
                            m_allow_changes[e * 2 + k] = 0;
                        }
                    }
                }
            }
        }

        build_integer_constraints();

        // Compute Max Flow
        hierarchy.DownsampleEdgeGraph(m_face_edge_orientation, m_face_edge_ids, m_edge_difference, m_allow_changes, 1);
        Optimizer::optimize_integer_constraints(hierarchy, m_singularities);
        hierarchy.UpdateGraphValue(m_face_edge_orientation, m_face_edge_ids, m_edge_difference);

        // potential bug
        subdivide_edge_to_length_considering_edge_differences(
            F, V, N, Q, O, &hierarchy.m_scales[0], m_V2E, hierarchy.m_E2E,
            m_boundary, m_non_manifold,
            m_edge_difference, m_edge_values, m_face_edge_orientation,
            m_face_edge_ids, m_edges_preserve,
            m_singularities, 1
        );

        m_allow_changes.clear();
        m_allow_changes.resize(m_edge_difference.size() * 2, 1);
        for (int i = 0; i < m_edges_preserve.size(); ++i) {
            if (m_edges_preserve[i] == 0) continue;
            int e = m_face_edge_ids[i / 3][i % 3];
            for (int k = 0; k < 2; ++k) {
                if (m_edge_difference[e][k] == 0) m_allow_changes[e * 2 + k] = 0;
            }
        }

        fix_flip_hierarchy();
        subdivide_edge_to_length_considering_edge_differences(
            F, V, N, Q, O, &hierarchy.m_scales[0], m_V2E, hierarchy.m_E2E,
            m_boundary, m_non_manifold,
            m_edge_difference, m_edge_values, m_face_edge_orientation,
            m_face_edge_ids, m_edges_preserve,
            m_singularities, 1
        );

        std::set<int> sharp_vertices;
        for (int i = 0; i < m_edges_preserve.size(); ++i) {
            if (m_edges_preserve[i] == 1) {
                sharp_vertices.insert(F(i % 3, i / 3));
                sharp_vertices.insert(F((i + 1) % 3, i / 3));
            }
        }

        Optimizer::optimize_positions_sharp(
            hierarchy,
            m_edge_values,
            m_edge_difference,
            m_edges_preserve,
            sharp_vertices,
            sharp_constraints,
            with_scale
        );

        Optimizer::optimize_positions_fixed(
            hierarchy,
            m_edge_values,
            m_edge_difference,
            sharp_vertices,
            sharp_constraints,
            with_scale
        );

        extract_quad();
        fix_valence();

        std::vector<int> sharp_o(m_positions_compact.size(), 0);
        std::map<int, std::pair<Vector3d, Vector3d> > compact_sharp_constraints;
        for (int i = 0; i < m_vertices_set.size(); ++i) {
            int sharpv = -1;
            for (auto &p: m_vertices_set[i]) {
                if (sharp_constraints.count(p)) {
                    sharpv = p;
                    sharp_o[i] = 1;
                    if (compact_sharp_constraints.count(i) == 0 ||
                        compact_sharp_constraints[i].second != Vector3d::Zero()) {
                        compact_sharp_constraints[i] = sharp_constraints[sharpv];
                        m_positions_compact[i] = O.col(sharpv);
                        compact_sharp_constraints[i].first = m_positions_compact[i];
                    }
                }
            }
        }

        std::map<std::pair<int, int>, int> o2e;
        for (int i = 0; i < m_faces_compact.size(); ++i) {
            for (int j = 0; j < 4; ++j) {
                int v1 = m_faces_compact[i][j];
                int v2 = m_faces_compact[i][(j + 1) % 4];
                o2e[std::make_pair(v1, v2)] = i * 4 + j;
            }
        }

        std::vector<std::vector<int> > v2o(V.cols());
        for (int i = 0; i < m_vertices_set.size(); ++i) {
            for (auto v: m_vertices_set[i]) {
                v2o[v].push_back(i);
            }
        }
        std::vector<Vector3d> diffs(m_faces_compact.size() * 4, Vector3d(0, 0, 0));
        std::vector<int> diff_count(m_faces_compact.size() * 4, 0);

        for (int i = 0; i < F.cols(); ++i) {
            for (int j = 0; j < 3; ++j) {
                int v1 = F(j, i);
                int v2 = F((j + 1) % 3, i);
                if (v1 != m_edge_values[m_face_edge_ids[i][j]].x) continue;
                if (m_edge_difference[m_face_edge_ids[i][j]].array().abs().sum() != 1) continue;
                if (v2o[v1].size() > 1 || v2o[v2].size() > 1) continue;
                for (auto o1: v2o[v1]) {
                    for (auto o2: v2o[v2]) {
                        auto key = std::make_pair(o1, o2);
                        if (o2e.count(key)) {
                            int dedge = o2e[key];
                            Vector3d q_1 = Q.col(v1);
                            Vector3d q_2 = Q.col(v2);
                            Vector3d n_1 = N.col(v1);
                            Vector3d n_2 = N.col(v2);
                            Vector3d q_1_y = n_1.cross(q_1);
                            Vector3d q_2_y = n_2.cross(q_2);
                            auto index = mathext::compat_orientation_extrinsic_index_4(q_1, n_1, q_2, n_2);
                            double s_x1 = S(0, v1), s_y1 = S(1, v1);
                            double s_x2 = S(0, v2), s_y2 = S(1, v2);
                            int rank_diff = (index.second + 4 - index.first) % 4;
                            if (rank_diff % 2 == 1) std::swap(s_x2, s_y2);
                            Vector3d qd_x = 0.5 * (mathext::rotate90_by(q_2, n_2, rank_diff) + q_1);
                            Vector3d qd_y = 0.5 * (mathext::rotate90_by(q_2_y, n_2, rank_diff) + q_1_y);
                            double scale_x = (with_scale ? 0.5 * (s_x1 + s_x2) : 1) * hierarchy.m_scale;
                            double scale_y = (with_scale ? 0.5 * (s_y1 + s_y2) : 1) * hierarchy.m_scale;
                            Vector2i diff = m_edge_difference[m_face_edge_ids[i][j]];
                            Vector3d C = diff[0] * scale_x * qd_x + diff[1] * scale_y * qd_y;

                            diff_count[dedge] += 1;
                            diffs[dedge] += C;
                            auto key = std::make_pair(o2, o1);
                            if (o2e.count(key)) {
                                int dedge = o2e[key];
                                diff_count[dedge] += 1;
                                diffs[dedge] -= C;
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < F.cols(); ++i) {
            Vector2i d1 = mathext::rshift90(m_edge_difference[m_face_edge_ids[i][0]], m_face_edge_orientation[i][0]);
            Vector2i d2 = mathext::rshift90(m_edge_difference[m_face_edge_ids[i][1]], m_face_edge_orientation[i][1]);
            if (d1[0] * d2[1] - d1[1] * d2[0] < 0) {
                for (int j = 0; j < 3; ++j) {
                    int v1 = F(j, i);
                    int v2 = F((j + 1) % 3, i);
                    for (auto o1: v2o[v1]) {
                        for (auto o2: v2o[v2]) {
                            auto key = std::make_pair(o1, o2);
                            if (o2e.count(key)) {
                                int dedge = o2e[key];
                                diff_count[dedge] = 0;
                                diffs[dedge] = Vector3d(0, 0, 0);
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < diff_count.size(); ++i) {
            if (diff_count[i] != 0) {
                diffs[i] /= diff_count[i];
                diff_count[i] = 1;
            }
        }

        Optimizer::optimize_positions_dynamic(
            F, V, N, Q,
            m_vertices_set, m_positions_compact, m_faces_compact, m_V2E_compact,
            m_E2E_compact,
            sqrt(m_surface_area / m_faces_compact.size()),
            diffs,
            diff_count,
            o2e,
            sharp_o,
            compact_sharp_constraints,
            with_scale
        );
    }

    // Post-Processing
    void Parametrizer::fix_valence() {
        // Remove Valence 2
        while (true) {
            bool update = false;
            std::vector<int> marks(m_V2E_compact.size(), 0);
            std::vector<int> erasedF(m_faces_compact.size(), 0);
            for (int i = 0; i < m_V2E_compact.size(); ++i) {
                int deid0 = m_V2E_compact[i];
                if (marks[i] || deid0 == -1) continue;
                int deid = deid0;
                std::vector<int> dedges;
                do {
                    dedges.push_back(deid);
                    int deid1 = deid / 4 * 4 + (deid + 3) % 4;
                    deid = m_E2E_compact[deid1];
                } while (deid != deid0 && deid != -1);
                if (dedges.size() == 2) {
                    int v1 = m_faces_compact[dedges[0] / 4][(dedges[0] + 1) % 4];
                    int v2 = m_faces_compact[dedges[0] / 4][(dedges[0] + 2) % 4];
                    int v3 = m_faces_compact[dedges[1] / 4][(dedges[1] + 1) % 4];
                    int v4 = m_faces_compact[dedges[1] / 4][(dedges[1] + 2) % 4];
                    if (marks[v1] || marks[v2] || marks[v3] || marks[v4]) continue;
                    marks[v1] = true;
                    marks[v2] = true;
                    marks[v3] = true;
                    marks[v4] = true;
                    if (v1 == v2 || v1 == v3 || v1 == v4 || v2 == v3 || v2 == v4 || v3 == v4) {
                        erasedF[dedges[0] / 4] = 1;
                    } else {
                        m_faces_compact[dedges[0] / 4] = Vector4i(v1, v2, v3, v4);
                    }
                    erasedF[dedges[1] / 4] = 1;
                    update = true;
                }
            }
            if (update) {
                int top = 0;
                for (int i = 0; i < erasedF.size(); ++i) {
                    if (erasedF[i] == 0) {
                        m_faces_compact[top++] = m_faces_compact[i];
                    }
                }
                m_faces_compact.resize(top);
                mathext::compute_direct_graph_quad(m_positions_compact, m_faces_compact, m_V2E_compact, m_E2E_compact,
                                                   m_boundary_compact, m_non_manifold_compact);
            } else {
                break;
            }
        }
        std::vector<std::vector<int> > v_dedges(m_V2E_compact.size());
        for (int i = 0; i < m_faces_compact.size(); ++i) {
            for (int j = 0; j < 4; ++j) {
                v_dedges[m_faces_compact[i][j]].push_back(i * 4 + j);
            }
        }
        int top = m_V2E_compact.size();
        for (int i = 0; i < v_dedges.size(); ++i) {
            std::map<int, int> groups;
            int group_id = 0;
            for (int j = 0; j < v_dedges[i].size(); ++j) {
                int deid = v_dedges[i][j];
                if (groups.count(deid)) continue;
                int deid0 = deid;
                do {
                    groups[deid] = group_id;
                    deid = deid / 4 * 4 + (deid + 3) % 4;
                    deid = m_E2E_compact[deid];
                } while (deid != deid0 && deid != -1);
                if (deid == -1) {
                    deid = deid0;
                    while (m_E2E_compact[deid] != -1) {
                        deid = m_E2E_compact[deid];
                        deid = deid / 4 * 4 + (deid + 1) % 4;
                        groups[deid] = group_id;
                    }
                }
                group_id += 1;
            }
            if (group_id > 1) {
                for (auto &g: groups) {
                    if (g.second >= 1) m_faces_compact[g.first / 4][g.first % 4] = top - 1 + g.second;
                }
                for (int j = 1; j < group_id; ++j) {
                    m_vertices_set.push_back(m_vertices_set[i]);
                    m_normals_compact.push_back(m_normals_compact[i]);
                    m_orientations_compact.push_back(m_orientations_compact[i]);
                    m_positions_compact.push_back(m_positions_compact[i]);
                }
                top = m_positions_compact.size();
            }
        }
        mathext::compute_direct_graph_quad(m_positions_compact, m_faces_compact, m_V2E_compact, m_E2E_compact,
                                           m_boundary_compact,
                                           m_non_manifold_compact);

        // Decrease Valence
        while (true) {
            bool update = false;
            std::vector<int> marks(m_V2E_compact.size(), 0);
            std::vector<int> valences(m_V2E_compact.size(), 0);
            for (int i = 0; i < m_V2E_compact.size(); ++i) {
                int deid0 = m_V2E_compact[i];
                if (deid0 == -1) continue;
                int deid = deid0;
                int count = 0;
                do {
                    count += 1;
                    int deid1 = m_E2E_compact[deid];
                    if (deid1 == -1) {
                        count += 1;
                        break;
                    }
                    deid = deid1 / 4 * 4 + (deid1 + 1) % 4;
                } while (deid != deid0 && deid != -1);
                if (deid == -1) count += 1;
                valences[i] = count;
            }
            std::priority_queue<std::pair<int, int> > prior_queue;
            for (int i = 0; i < valences.size(); ++i) {
                if (valences[i] > 5) prior_queue.push(std::make_pair(valences[i], i));
            }
            while (!prior_queue.empty()) {
                auto info = prior_queue.top();
                prior_queue.pop();
                if (marks[info.second]) continue;
                int deid0 = m_V2E_compact[info.second];
                if (deid0 == -1) continue;
                int deid = deid0;
                std::vector<int> loop_vertices, loop_dedges;;
                bool marked = false;
                do {
                    int v = m_faces_compact[deid / 4][(deid + 1) % 4];
                    loop_dedges.push_back(deid);
                    loop_vertices.push_back(v);
                    if (marks[v]) marked = true;
                    int deid1 = m_E2E_compact[deid];
                    if (deid1 == -1) break;
                    deid = deid1 / 4 * 4 + (deid1 + 1) % 4;
                } while (deid != deid0 && deid != -1);
                if (marked) continue;

                if (deid != -1) {
                    int step = (info.first + 1) / 2;
                    std::pair<int, int> min_val(0x7fffffff, 0x7fffffff);
                    int split_idx = -1;
                    for (int i = 0; i < loop_vertices.size(); ++i) {
                        if (i + step >= loop_vertices.size()) continue;
                        int v1 = valences[loop_vertices[i]];
                        int v2 = valences[loop_vertices[i + step]];
                        if (v1 < v2) std::swap(v1, v2);
                        auto key = std::make_pair(v1, v2);
                        if (key < min_val) {
                            min_val = key;
                            split_idx = i + 1;
                        }
                    }
                    if (min_val.first >= info.first) continue;
                    update = true;
                    for (int id = split_idx; id < split_idx + step; ++id) {
                        m_faces_compact[loop_dedges[id] / 4][loop_dedges[id] % 4] = m_positions_compact.size();
                    }
                    m_faces_compact.push_back(Vector4i(
                        m_positions_compact.size(),
                        loop_vertices[(split_idx + loop_vertices.size() - 1) % loop_vertices.size()],
                        info.second,
                        loop_vertices[(split_idx + step - 1 + loop_vertices.size()) %
                                      loop_vertices.size()]));
                } else {
                    for (int id = loop_vertices.size() / 2; id < loop_vertices.size(); ++id) {
                        m_faces_compact[loop_dedges[id] / 4][loop_dedges[id] % 4] = m_positions_compact.size();
                    }
                    update = true;
                }
                marks[info.second] = 1;
                for (int i = 0; i < loop_vertices.size(); ++i) {
                    marks[loop_vertices[i]] = 1;
                }
                m_vertices_set.push_back(m_vertices_set[info.second]);
                m_positions_compact.push_back(m_positions_compact[info.second]);
                m_normals_compact.push_back(m_normals_compact[info.second]);
                m_orientations_compact.push_back(m_orientations_compact[info.second]);
            }
            if (!update) {
                break;
            } else {
                mathext::compute_direct_graph_quad(m_positions_compact, m_faces_compact, m_V2E_compact, m_E2E_compact,
                                                   m_boundary_compact, m_non_manifold_compact);
            }
        }

        // Remove Zero Valence
        std::vector<int> valences(m_V2E_compact.size(), 0);
        for (int i = 0; i < m_faces_compact.size(); ++i) {
            for (int j = 0; j < 4; ++j) {
                valences[m_faces_compact[i][j]] = 1;
            }
        }
        top = 0;
        std::vector<int> compact_indices(valences.size());
        for (int i = 0; i < valences.size(); ++i) {
            if (valences[i] == 0) continue;
            m_normals_compact[top] = m_normals_compact[i];
            m_positions_compact[top] = m_positions_compact[i];
            m_orientations_compact[top] = m_orientations_compact[i];
            m_vertices_set[top] = m_vertices_set[i];
            compact_indices[i] = top;
            top += 1;
        }
        for (int i = 0; i < m_faces_compact.size(); ++i) {
            for (int j = 0; j < 4; ++j) {
                m_faces_compact[i][j] = compact_indices[m_faces_compact[i][j]];
            }
        }
        m_normals_compact.resize(top);
        m_positions_compact.resize(top);
        m_orientations_compact.resize(top);
        m_vertices_set.resize(top);
        mathext::compute_direct_graph_quad(m_positions_compact, m_faces_compact, m_V2E_compact, m_E2E_compact,
                                           m_boundary_compact,
                                           m_non_manifold_compact); {
            mathext::compute_direct_graph_quad(m_positions_compact, m_faces_compact, m_V2E_compact, m_E2E_compact,
                                               m_boundary_compact,
                                               m_non_manifold_compact);
            std::vector<int> masks(m_faces_compact.size() * 4, 0);
            for (int i = 0; i < m_V2E_compact.size(); ++i) {
                int deid0 = m_V2E_compact[i];
                if (deid0 == -1) continue;
                int deid = deid0;
                do {
                    masks[deid] = 1;
                    deid = m_E2E_compact[deid];
                    if (deid == -1) {
                        break;
                    }
                    deid = deid / 4 * 4 + (deid + 1) % 4;
                } while (deid != deid0 && deid != -1);
            }
            std::vector<std::vector<int> > v_dedges(m_V2E_compact.size());
            for (int i = 0; i < m_faces_compact.size(); ++i) {
                for (int j = 0; j < 4; ++j) {
                    v_dedges[m_faces_compact[i][j]].push_back(i * 4 + j);
                }
            }
        }
        std::map<int, int> pts;
        for (int i = 0; i < m_V2E_compact.size(); ++i) {
            int deid0 = m_V2E_compact[i];
            if (deid0 == -1) continue;
            int deid = deid0;
            int count = 0;
            do {
                count += 1;
                int deid1 = m_E2E_compact[deid];
                if (deid1 == -1) break;
                deid = deid1 / 4 * 4 + (deid1 + 1) % 4;
            } while (deid != deid0 && deid != -1);
            if (pts.count(count) == 0)
                pts[count] = 1;
            else
                pts[count] += 1;
        }
    }

    void Parametrizer::fix_flip_hierarchy() {
        Hierarchy fh;
        fh.DownsampleEdgeGraph(m_face_edge_orientation, m_face_edge_ids, m_edge_difference, m_allow_changes, -1);
        fh.FixFlip();
        fh.UpdateGraphValue(m_face_edge_orientation, m_face_edge_ids, m_edge_difference);
    }

    void Parametrizer::find_fix_holes() {
        for (int i = 0; i < m_faces_compact.size(); ++i) {
            for (int j = 0; j < 4; ++j) {
                int v1 = m_faces_compact[i][j];
                int v2 = m_faces_compact[i][(j + 1) % 4];
                auto key = std::make_pair(v1, v2);
                m_quad_edges.insert(key);
            }
        }
        std::vector<int> detected_boundary(m_E2E_compact.size(), 0);
        for (int i = 0; i < m_E2E_compact.size(); ++i) {
            if (detected_boundary[i] != 0 || m_E2E_compact[i] != -1) continue;
            std::vector<int> loop_edges;
            int current_e = i;

            while (detected_boundary[current_e] == 0) {
                detected_boundary[current_e] = 1;
                loop_edges.push_back(current_e);
                current_e = current_e / 4 * 4 + (current_e + 1) % 4;
                while (m_E2E_compact[current_e] != -1) {
                    current_e = m_E2E_compact[current_e];
                    current_e = current_e / 4 * 4 + (current_e + 1) % 4;
                }
            }
            std::vector<int> loop_vertices(loop_edges.size());
            for (int j = 0; j < loop_edges.size(); ++j) {
                loop_vertices[j] = m_faces_compact[loop_edges[j] / 4][loop_edges[j] % 4];
            }
            if (loop_vertices.size() < 25) close_hole(loop_vertices);
        }
    }

    void Parametrizer::close_hole(std::vector<int> &loop_vertices) {
        std::vector<std::vector<int> > loop_vertices_array;
        std::unordered_map<int, int> map_loops;
        for (int i = 0; i < loop_vertices.size(); ++i) {
            if (map_loops.count(loop_vertices[i])) {
                int j = map_loops[loop_vertices[i]];
                loop_vertices_array.push_back(std::vector<int>());
                if (i - j > 3 && (i - j) % 2 == 0) {
                    for (int k = j; k < i; ++k) {
                        if (map_loops.count(loop_vertices[k])) {
                            loop_vertices_array.back().push_back(loop_vertices[k]);
                            map_loops.erase(loop_vertices[k]);
                        }
                    }
                }
            }
            map_loops[loop_vertices[i]] = i;
        }
        if (map_loops.size() >= 3) {
            loop_vertices_array.push_back(std::vector<int>());
            for (int k = 0; k < loop_vertices.size(); ++k) {
                if (map_loops.count(loop_vertices[k])) {
                    if (map_loops.count(loop_vertices[k])) {
                        loop_vertices_array.back().push_back(loop_vertices[k]);
                        map_loops.erase(loop_vertices[k]);
                    }
                }
            }
        }
        for (int i = 0; i < loop_vertices_array.size(); ++i) {
            auto &loop_vertices = loop_vertices_array[i];
            if (loop_vertices.size() == 0) return;
            std::vector<Vector4i> quads;
            compute_quad_energy(loop_vertices, quads, 0);
            for (auto &p: quads) {
                bool flag = false;
                for (int j = 0; j < 4; ++j) {
                    int v1 = p[j];
                    int v2 = p[(j + 1) % 4];
                    auto key = std::make_pair(v1, v2);
                    if (m_quad_edges.count(key)) {
                        flag = true;
                        break;
                    }
                }
                if (!flag) {
                    for (int j = 0; j < 4; ++j) {
                        int v1 = p[j];
                        int v2 = p[(j + 1) % 4];
                        auto key = std::make_pair(v1, v2);
                        m_quad_edges.insert(key);
                    }
                    m_faces_compact.push_back(p);
                }
            }
        }
    }

    double Parametrizer::compute_quad_energy(
        std::vector<int> &loop_vertices,
        std::vector<Vector4i> &res_quads,
        int level
    ) {
        if (loop_vertices.size() < 4) return 0;

        if (loop_vertices.size() == 4) {
            double energy = 0;
            for (int j = 0; j < 4; ++j) {
                int v0 = loop_vertices[j];
                int v2 = loop_vertices[(j + 1) % 4];
                int v1 = loop_vertices[(j + 3) % 4];
                Vector3d pt1 = (m_positions_compact[v1] - m_positions_compact[v0]).normalized();
                Vector3d pt2 = (m_positions_compact[v2] - m_positions_compact[v0]).normalized();
                Vector3d n = pt1.cross(pt2);
                double sina = n.norm();
                if (n.dot(m_normals_compact[v0]) < 0) sina = -sina;
                double cosa = pt1.dot(pt2);
                double angle = atan2(sina, cosa) / 3.141592654 * 180.0;
                if (angle < 0) angle = 360 + angle;
                energy += angle * angle;
            }
            res_quads.push_back(
                Vector4i(loop_vertices[0], loop_vertices[3], loop_vertices[2], loop_vertices[1]));
            return energy;
        }
        double max_energy = 1e30;
        for (int seg1 = 2; seg1 < loop_vertices.size(); seg1 += 2) {
            for (int seg2 = seg1 + 1; seg2 < loop_vertices.size(); seg2 += 2) {
                std::vector<Vector4i> quads[4];
                std::vector<int> vertices = {
                    loop_vertices[0], loop_vertices[1], loop_vertices[seg1],
                    loop_vertices[seg2]
                };
                double energy = 0;
                energy += compute_quad_energy(vertices, quads[0], level + 1);
                if (seg1 > 2) {
                    std::vector<int> vertices(loop_vertices.begin() + 1, loop_vertices.begin() + seg1);
                    vertices.push_back(loop_vertices[seg1]);
                    energy += compute_quad_energy(vertices, quads[1], level + 1);
                }
                if (seg2 != seg1 + 1) {
                    std::vector<int> vertices(loop_vertices.begin() + seg1,
                                              loop_vertices.begin() + seg2);
                    vertices.push_back(loop_vertices[seg2]);
                    energy += compute_quad_energy(vertices, quads[2], level + 2);
                }
                if (seg2 + 1 != loop_vertices.size()) {
                    std::vector<int> vertices(loop_vertices.begin() + seg2, loop_vertices.end());
                    vertices.push_back(loop_vertices[0]);
                    energy += compute_quad_energy(vertices, quads[3], level + 1);
                }
                if (max_energy > energy) {
                    max_energy = energy;
                    res_quads.clear();
                    for (int i = 0; i < 4; ++i) {
                        for (auto &v: quads[i]) {
                            res_quads.push_back(v);
                        }
                    }
                }
            }
        }
        return max_energy;
    }

    void Parametrizer::extract_quad() {
        Hierarchy fh;
        fh.DownsampleEdgeGraph(m_face_edge_orientation, m_face_edge_ids, m_edge_difference, m_allow_changes, -1);
        auto &V = m_hierarchy.m_vertices[0];
        auto &F = m_hierarchy.m_faces;
        m_disjoint_tree = DisjointTree(V.cols());
        auto &diffs = fh.mEdgeDiff.front();
        for (int i = 0; i < diffs.size(); ++i) {
            if (diffs[i] == Vector2i::Zero()) {
                m_disjoint_tree.Merge(m_edge_values[i].x, m_edge_values[i].y);
            }
        }
        m_disjoint_tree.BuildCompactParent();
        auto &F2E = fh.mF2E.back();
        auto &E2F = fh.mE2F.back();
        auto &EdgeDiff = fh.mEdgeDiff.back();
        auto &FQ = fh.mFQ.back();

        std::vector<int> edge(E2F.size());
        std::vector<int> face(F2E.size());
        for (int i = 0; i < diffs.size(); ++i) {
            int t = i;
            for (int j = 0; j < fh.mToUpperEdges.size(); ++j) {
                t = fh.mToUpperEdges[j][t];
                if (t < 0) break;
            }
            if (t >= 0) edge[t] = i;
        }
        for (int i = 0; i < F.cols(); ++i) {
            int t = i;
            for (int j = 0; j < fh.mToUpperFaces.size(); ++j) {
                t = fh.mToUpperFaces[j][t];
                if (t < 0) break;
            }
            if (t >= 0) face[t] = i;
        }
        fh.UpdateGraphValue(m_face_edge_orientation, m_face_edge_ids, m_edge_difference);

        auto &O = m_hierarchy.m_positions[0];
        auto &Q = m_hierarchy.m_orientation[0];
        auto &N = m_hierarchy.m_normals[0];
        int num_v = m_disjoint_tree.CompactNum();
        m_vertices_set.resize(num_v);
        m_positions_compact.resize(num_v, Vector3d::Zero());
        m_orientations_compact.resize(num_v, Vector3d::Zero());
        m_normals_compact.resize(num_v, Vector3d::Zero());
        m_counter.resize(num_v, 0);
        for (int i = 0; i < O.cols(); ++i) {
            int compact_v = m_disjoint_tree.Index(i);
            m_vertices_set[compact_v].push_back(i);
            m_positions_compact[compact_v] += O.col(i);
            m_normals_compact[compact_v] = m_normals_compact[compact_v] * m_counter[compact_v] + N.col(i);
            m_normals_compact[compact_v].normalize();
            if (m_counter[compact_v] == 0)
                m_orientations_compact[compact_v] = Q.col(i);
            else {
                auto pairs = mathext::compat_orientation_extrinsic_4(m_orientations_compact[compact_v],
                                                                     m_normals_compact[compact_v],
                                                                     Q.col(i), N.col(i));
                m_orientations_compact[compact_v] = (pairs.first * m_counter[compact_v] + pairs.second).normalized();
            }
            m_counter[compact_v] += 1;
        }
        for (int i = 0; i < m_positions_compact.size(); ++i) {
            m_positions_compact[i] /= m_counter[i];
        }

        build_triangle_manifold(m_disjoint_tree, edge, face, m_edge_values, F2E, E2F, EdgeDiff, FQ);
    }

    void Parametrizer::build_triangle_manifold(
        DisjointTree &disajoint_tree,
        std::vector<int> &edge,
        std::vector<int> &face,
        std::vector<DEdge> &edge_values,
        std::vector<Vector3i> &F2E,
        std::vector<Vector2i> &E2F,
        std::vector<Vector2i> &EdgeDiff,
        std::vector<Vector3i> &FQ
    ) {
        auto &F = m_hierarchy.m_faces;
        std::vector<int> E2E(F2E.size() * 3, -1);
        for (int i = 0; i < E2F.size(); ++i) {
            int v1 = E2F[i][0];
            int v2 = E2F[i][1];
            int t1 = 0;
            int t2 = 2;
            if (v1 != -1)
                while (F2E[v1][t1] != i) t1 += 1;
            if (v2 != -1)
                while (F2E[v2][t2] != i) t2 -= 1;
            t1 += v1 * 3;
            t2 += v2 * 3;
            if (v1 != -1)
                E2E[t1] = (v2 == -1) ? -1 : t2;
            if (v2 != -1)
                E2E[t2] = (v1 == -1) ? -1 : t1;
        }

        std::vector<Vector3i> triangle_vertices(F2E.size(), Vector3i(-1, -1, -1));
        int num_v = 0;
        std::vector<Vector3d> N, Q, O;
        std::vector<std::vector<int> > Vs;
        for (int i = 0; i < F2E.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                if (triangle_vertices[i][j] != -1) continue;
                int f = face[i];
                int v = disajoint_tree.Index(F(j, f));
                Vs.push_back(m_vertices_set[v]);
                Q.push_back(m_orientations_compact[v]);
                N.push_back(m_normals_compact[v]);
                O.push_back(m_positions_compact[v]);
                int deid0 = i * 3 + j;
                int deid = deid0;
                do {
                    triangle_vertices[deid / 3][deid % 3] = num_v;
                    deid = E2E[deid / 3 * 3 + (deid + 2) % 3];
                } while (deid != deid0 && deid != -1);
                if (deid == -1) {
                    deid = deid0;
                    do {
                        deid = E2E[deid];
                        if (deid == -1)
                            break;
                        deid = deid / 3 * 3 + (deid + 1) % 3;
                        triangle_vertices[deid / 3][deid % 3] = num_v;
                    } while (deid != -1);
                }
                num_v += 1;
            }
        }

        int num_v0 = num_v;
        do {
            num_v0 = num_v;
            std::vector<std::vector<int> > vert_to_dedge(num_v);
            for (int i = 0; i < triangle_vertices.size(); ++i) {
                Vector3i pt = triangle_vertices[i];
                if (pt[0] == pt[1] || pt[1] == pt[2] || pt[2] == pt[0]) {
                    for (int j = 0; j < 3; ++j) {
                        int t = E2E[i * 3 + j];
                        if (t != -1) E2E[t] = -1;
                    }
                    for (int j = 0; j < 3; ++j) {
                        E2E[i * 3 + j] = -1;
                    }
                } else {
                    for (int j = 0; j < 3; ++j)
                        vert_to_dedge[triangle_vertices[i][j]].push_back(i * 3 + j);
                }
            }
            std::vector<int> colors(triangle_vertices.size() * 3, -1),
                    reverse_colors(triangle_vertices.size() * 3, -1);
            for (int i = 0; i < vert_to_dedge.size(); ++i) {
                int num_color = 0;
                for (int j = 0; j < vert_to_dedge[i].size(); ++j) {
                    int deid = vert_to_dedge[i][j];
                    if (colors[deid] != -1) continue;
                    std::list<int> l;
                    int deid0 = deid;
                    do {
                        l.push_back(deid);
                        deid = deid / 3 * 3 + (deid + 2) % 3;
                        deid = E2E[deid];
                    } while (deid != -1 && deid != deid0);
                    if (deid == -1) {
                        deid = deid0;
                        do {
                            deid = E2E[deid];
                            if (deid == -1) break;
                            deid = deid / 3 * 3 + (deid + 1) % 3;
                            if (deid == deid0) break;
                            l.push_front(deid);
                        } while (true);
                    }
                    std::vector<int> dedges;
                    for (auto &e: l) dedges.push_back(e);
                    std::map<std::pair<int, int>, int> loc;
                    std::vector<int> deid_colors(dedges.size(), num_color);
                    num_color += 1;
                    for (int jj = 0; jj < dedges.size(); ++jj) {
                        int deid = dedges[jj];
                        colors[deid] = 0;
                        int v1 = triangle_vertices[deid / 3][deid % 3];
                        int v2 = triangle_vertices[deid / 3][(deid + 1) % 3];
                        std::pair<int, int> pt(v1, v2);
                        if (loc.count(pt)) {
                            int s = loc[pt];
                            for (int k = s; k < jj; ++k) {
                                int deid1 = dedges[k];
                                int v11 = triangle_vertices[deid1 / 3][deid1 % 3];
                                int v12 = triangle_vertices[deid1 / 3][(deid1 + 1) % 3];
                                std::pair<int, int> pt1(v11, v12);
                                loc.erase(pt1);
                                deid_colors[k] = num_color;
                            }
                            num_color += 1;
                        }
                        loc[pt] = jj;
                    }
                    for (int j = 0; j < dedges.size(); ++j) {
                        int deid = dedges[j];
                        int color = deid_colors[j];
                        if (color > 0) {
                            triangle_vertices[deid / 3][deid % 3] = num_v + color - 1;
                        }
                    }
                }
                if (num_color > 1) {
                    for (int j = 0; j < num_color - 1; ++j) {
                        Vs.push_back(Vs[i]);
                        O.push_back(O[i]);
                        N.push_back(N[i]);
                        Q.push_back(Q[i]);
                    }
                    num_v += num_color - 1;
                }
            }
        } while (num_v != num_v0);
        int offset = 0;
        std::vector<Vector3i> triangle_edges, triangle_orients;
        for (int i = 0; i < triangle_vertices.size(); ++i) {
            Vector3i pt = triangle_vertices[i];
            if (pt[0] == pt[1] || pt[1] == pt[2] || pt[2] == pt[0]) continue;
            triangle_vertices[offset++] = triangle_vertices[i];
            triangle_edges.push_back(F2E[i]);
            triangle_orients.push_back(FQ[i]);
        }
        triangle_vertices.resize(offset);
        std::set<int> flip_vertices;
        for (int i = 0; i < triangle_vertices.size(); ++i) {
            Vector2i d1 = mathext::rshift90(EdgeDiff[triangle_edges[i][0]], triangle_orients[i][0]);
            Vector2i d2 = mathext::rshift90(EdgeDiff[triangle_edges[i][1]], triangle_orients[i][1]);
            int area = d1[0] * d2[1] - d1[1] * d2[0];
            if (area < 0) {
                for (int j = 0; j < 3; ++j) {
                    flip_vertices.insert(triangle_vertices[i][j]);
                }
            }
        }
        MatrixXd NV(3, num_v);
        MatrixXi NF(3, triangle_vertices.size());
        memcpy(NF.data(), triangle_vertices.data(), sizeof(int) * 3 * triangle_vertices.size());
        VectorXi NV2E, NE2E, NB, NN;
        mathext::compute_direct_graph(NV, NF, NV2E, NE2E, NB, NN);

        std::map<DEdge, std::pair<Vector3i, Vector3i> > quads;
        for (int i = 0; i < triangle_vertices.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                int e = triangle_edges[i][j];
                int v1 = triangle_vertices[i][j];
                int v2 = triangle_vertices[i][(j + 1) % 3];
                int v3 = triangle_vertices[i][(j + 2) % 3];
                if (abs(EdgeDiff[e][0]) == 1 && abs(EdgeDiff[e][1]) == 1) {
                    DEdge edge(v1, v2);
                    if (quads.count(edge))
                        quads[edge].second = Vector3i(v1, v2, v3);
                    else
                        quads[edge] = std::make_pair(Vector3i(v1, v2, v3), Vector3i(-1, -1, -1));
                }
            }
        }

        for (auto &p: quads) {
            if (p.second.second[0] != -1 && p.second.first[2] != p.second.second[2]) {
                m_faces_compact.push_back(Vector4i(p.second.first[1], p.second.first[2], p.second.first[0],
                                                   p.second.second[2]));
            }
        }
        std::swap(Vs, m_vertices_set);
        std::swap(m_positions_compact, O);
        std::swap(m_normals_compact, N);
        std::swap(m_orientations_compact, Q);

        mathext::compute_direct_graph_quad(
            m_positions_compact,
            m_faces_compact,
            m_V2E_compact,
            m_E2E_compact,
            m_boundary_compact,
            m_non_manifold_compact
        );

        while (true) {
            std::vector<int> erasedF(m_faces_compact.size(), 0);
            for (int i = 0; i < m_faces_compact.size(); ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = j + 1; k < 4; ++k) {
                        if (m_faces_compact[i][j] == m_faces_compact[i][k]) {
                            erasedF[i] = 1;
                        }
                    }
                }
            }
            for (int i = 0; i < m_positions_compact.size(); ++i) {
                int v = 0;
                int e0 = m_V2E_compact[i];
                if (e0 == -1) continue;
                std::vector<int> dedges;
                int e = e0;
                do {
                    dedges.push_back(e);
                    v += 1;
                    e = e / 4 * 4 + (e + 3) % 4;
                    e = m_E2E_compact[e];
                } while (e != e0 && e != -1);
                if (e == -1) {
                    int e = e0;
                    while (true) {
                        e = m_E2E_compact[e];
                        if (e == -1) break;
                        e = e / 4 * 4 + (e + 1) % 4;
                        v += 1;
                        dedges.push_back(e);
                    }
                }
                if (v == 2) {
                    //                erasedF[dedges[1] / 4] = 1;
                    //                m_faces_compact[dedges[0]/4][dedges[0]%4] =
                    //                m_faces_compact[dedges[1]/4][(dedges[1]+2)%4];
                }
            }
            offset = 0;
            for (int i = 0; i < m_faces_compact.size(); ++i) {
                if (erasedF[i] == 0) m_faces_compact[offset++] = m_faces_compact[i];
            }
            if (offset == m_faces_compact.size()) break;
            m_faces_compact.resize(offset);
            mathext::compute_direct_graph_quad(
                m_positions_compact,
                m_faces_compact,
                m_V2E_compact,
                m_E2E_compact,
                m_boundary_compact,
                m_non_manifold_compact
            );
        }
        find_fix_holes();
        mathext::compute_direct_graph_quad(
            m_positions_compact,
            m_faces_compact,
            m_V2E_compact,
            m_E2E_compact,
            m_boundary_compact,
            m_non_manifold_compact
        );
    }
}
