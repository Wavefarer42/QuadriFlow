#include <algorithm>
#include <queue>
#include <unordered_map>

#include "hierarchy.h"
#include "field-math.h"
#include "pcg32.h"

namespace services {

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
            DownsampleGraph(m_adjacency[i], m_vertices[i], m_normals[i], m_vertex_area[i], m_vertices[i + 1], m_normals[i + 1], m_vertex_area[i + 1], mToUpper[i],
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
                coordinate_system(m_normals[i].col(j), s, t);
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
            const entities::AdjacentMatrix &adj,
            int size,
            std::vector<std::vector<int>> &phases
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
    ) {
        struct Entry {
            int i, j;
            double order;

            inline Entry() { i = j = -1; };

            inline Entry(int i, int j, double order) : i(i), j(j), order(order) {}

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
        std::vector<std::vector<entities::Link>> scratches(V_p.cols());
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
                for (auto &link: ad) scratch.push_back(entities::Link(to_lower[link.id], link.weight));
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
                diff += rshift90(edge_diff[F2E[i][j]], FQ[i][j]);
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
                    std::vector<std::pair<int, int>> paths;
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
            Vector2i diff1 = rshift90(EdgeDiff[F2E[f][0]], FQ[f][0]);
            Vector2i diff2 = rshift90(EdgeDiff[F2E[f][1]], FQ[f][1]);
            return diff1[0] * diff2[1] - diff1[1] * diff2[0];
        };
        std::vector<int> valences(F2E.size() * 3, -10000);  // comment this line
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
                diff = -rshift90(diff, FQ[deid / 3][deid % 3]);
                deid = deid / 3 * 3 + (deid + 1) % 3;
                // transform to local
                diff = rshift90(diff, (4 - FQ[deid / 3][deid % 3]) % 4);
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
                    nEdgeDiff[i] = rshift90(EdgeDiff[toUpper[i]], orient);
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
                    auto result = compat_orientation_extrinsic_4(q_i, n_i, q_j, n_j);
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
                    auto result = compat_position_extrinsic_4(
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
