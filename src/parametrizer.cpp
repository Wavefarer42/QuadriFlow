#include <list>
#include <map>
#include <fstream>
#include <set>
#include <random>

#include <Eigen/Sparse>

#include "services.h"
#include "field-math.h"
#include "optimizer.h"
#include "dedge.h"


namespace services {

    void Parametrizer::save_to_obj(const char *obj_name) {
        std::ofstream os(obj_name);
        for (int i = 0; i < m_positions_compact.size(); ++i) {
            auto t = m_positions_compact[i] * this->m_normalize_scale + this->m_normalize_offset;
            os << "v " << t[0] << " " << t[1] << " " << t[2] << "\n";
        }
        for (int i = 0; i < m_faces_compact.size(); ++i) {
            os << "f " << m_faces_compact[i][0] + 1 << " " << m_faces_compact[i][1] + 1 << " "
               << m_faces_compact[i][2] + 1 << " " << m_faces_compact[i][3] + 1 << "\n";
        }
        os.close();
    }

    void Parametrizer::normalize_mesh() {
        double maxV[3] = {-1e30, -1e30, -1e30};
        double minV[3] = {1e30, 1e30, 1e30};

        for (int i = 0; i < m_vertices.cols(); ++i) {
            for (int j = 0; j < 3; ++j) {
                maxV[j] = std::max(maxV[j], m_vertices(j, i));
                minV[j] = std::min(minV[j], m_vertices(j, i));
            }
        }
        double scale = std::max(std::max(maxV[0] - minV[0], maxV[1] - minV[1]), maxV[2] - minV[2]) * 0.5;
        for (int i = 0; i < m_vertices.cols(); ++i) {
            for (int j = 0; j < 3; ++j) {
                m_vertices(j, i) = (m_vertices(j, i) - (maxV[j] + minV[j]) * 0.5) / scale;
            }
        }
        this->m_normalize_scale = scale;
        this->m_normalize_offset = Vector3d(0.5 * (maxV[0] + minV[0]),
                                          0.5 * (maxV[1] + minV[1]),
                                          0.5 * (maxV[2] + minV[2]));
    }

    void Parametrizer::analyze_mesh() {
        m_surface_area = 0;
        m_average_edge_length = 0;
        m_max_edge_length = 0;
        for (int f = 0; f < m_faces.cols(); ++f) {
            Vector3d v[3] = {m_vertices.col(m_faces(0, f)), m_vertices.col(m_faces(1, f)),
                             m_vertices.col(m_faces(2, f))};
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
                int ep = dedge_prev_3(edge), en = dedge_next_3(edge);

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
                edge = dedge_next_3(opp);
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
                if (m_sharp_edges[edge]) break;
                edge = m_E2E[edge];
                if (edge != -1) edge = dedge_next_3(edge);
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
                double angle = fast_acos(d0.dot(d1) / std::sqrt(d0.squaredNorm() * d1.squaredNorm()));

                /* "Computing Vertex Normals from Polygonal Facets"
                 by Grit Thuermer and Charles A. Wuethrich, JGT 1998, Vol 3 */
                if (std::isfinite(angle)) normal += m_faces_normals.col(edge / 3) * angle;

                int opp = m_E2E[edge];
                if (opp == -1) break;

                edge = dedge_next_3(opp);
                if (m_sharp_edges[edge]) break;
            } while (edge != stop);
            double norm = normal.norm();
            m_normals_vertices.col(i) = norm > RCPOVERFLOW ? Vector3d(normal / norm) : Vector3d::UnitX();
        }
    }

    void Parametrizer::find_edges_and_features_and_boundaries(
            bool should_preserve_boundaries,
            bool should_preserve_edges
    ) {
        m_sharp_edges.resize(m_faces.cols() * 3, 0);

        if (should_preserve_boundaries) {
            for (int i = 0; i < m_sharp_edges.size(); ++i) {
                int re = m_E2E[i];
                if (re == -1) {
                    m_sharp_edges[i] = 1;
                }
            }
        }

        if (should_preserve_edges) return;

        std::vector<Vector3d> face_normals(m_faces.cols());
        for (int i = 0; i < m_faces.cols(); ++i) {
            Vector3d p1 = m_vertices.col(m_faces(0, i));
            Vector3d p2 = m_vertices.col(m_faces(1, i));
            Vector3d p3 = m_vertices.col(m_faces(2, i));
            face_normals[i] = (p2 - p1).cross(p3 - p1).normalized();
        }

        double cos_thres = cos(60.0 / 180.0 * 3.141592654);
        for (int i = 0; i < m_sharp_edges.size(); ++i) {
            int e = i;
            int re = m_E2E[e];
            Vector3d &n1 = face_normals[e / 3];
            Vector3d &n2 = face_normals[re / 3];
            if (n1.dot(n2) < cos_thres) {
                m_sharp_edges[i] = 1;
            }
        }
    }

    void generate_adjacency_matrix_uniform(
            const MatrixXi &F,
            const VectorXi &V2E,
            const VectorXi &E2E,
            const VectorXi &nonManifold,
            entities::AdjacentMatrix &adj
    ) {
        adj.resize(V2E.size());
        for (int i = 0; i < adj.size(); ++i) {
            int start = V2E[i];
            int edge = start;
            if (start == -1)
                continue;
            do {
                int base = edge % 3, f = edge / 3;
                int opp = E2E[edge], next = dedge_next_3(opp);
                if (adj[i].empty())
                    adj[i].push_back(entities::Link(F((base + 2) % 3, f)));
                if (opp == -1 || next != start) {
                    adj[i].push_back(entities::Link(F((base + 1) % 3, f)));
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
            const int e0p = E2E[dedge_prev_3(e0)], e0n = E2E[dedge_next_3(e0)];

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
                const int e1p = E2E[dedge_prev_3(e1)], e1n = E2E[dedge_next_3(e1)];
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

    void Parametrizer::initialize_parameterizer(
            bool should_preserve_boundaries,
            bool should_preserve_edges,
            int target_face_count,
            bool with_scale
    ) {
        m_hierarchy.clearConstraints();
        normalize_mesh();
        analyze_mesh();

        // initialize m_rho
        m_rho.resize(m_vertices.cols(), 1);
        for (int i = 0; i < m_vertices.cols(); ++i) {
            m_rho[i] = 1;
        }

        // initialize the scale of the mesh
        if (target_face_count <= 0) {
            m_scale = sqrt(m_surface_area / m_vertices.cols());
        } else {
            m_scale = std::sqrt(m_surface_area / target_face_count);
        }

        // Computes the directed graph and subdivides if the scale is larger than the maximum edge length.
        double target_len = std::min(m_scale / 2, m_average_edge_length * 2);
        if (target_len < m_max_edge_length) {
            while (!compute_direct_graph(m_vertices, m_faces, m_V2E, m_E2E, m_boundary, m_non_manifold));
            subdivide_edges_to_length(m_faces, m_vertices, m_rho, m_V2E, m_E2E, m_boundary, m_non_manifold, target_len);
        }
        while (!compute_direct_graph(m_vertices, m_faces, m_V2E, m_E2E, m_boundary, m_non_manifold));

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

        find_edges_and_features_and_boundaries(should_preserve_edges, should_preserve_boundaries);

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
                entities::DEdge e2(v1, v2);
                Vector2i diff2;
                int rank2;
                if (v1 > v2) {
                    rank2 = m_singularity_rank(k2, i);
                    diff2 = rshift90(
                            Vector2i(
                                    -m_singularity_index(k1 * 2, i),
                                    -m_singularity_index(k1 * 2 + 1, i)
                            ),
                            rank2
                    );
                } else {
                    rank2 = m_singularity_rank(k1, i);
                    diff2 = rshift90(
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
        std::vector<std::pair<int, int>> E2D(m_edge_difference.size(), std::make_pair(-1, -1));
        for (int i = 0; i < F.cols(); ++i) {
            int v0 = F(0, i);
            int v1 = F(1, i);
            int v2 = F(2, i);
            entities::DEdge e0(v0, v1), e1(v1, v2), e2(v2, v0);
            const Vector3i &eid = m_face_edge_ids[i];
            Vector2i variable_id[3];
            for (int i = 0; i < 3; ++i) {
                variable_id[i] = Vector2i(eid[i] * 2 + 1, eid[i] * 2 + 2);
            }
            auto index1 =
                    compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v1), N.col(v1));
            auto index2 =
                    compat_orientation_extrinsic_index_4(Q.col(v0), N.col(v0), Q.col(v2), N.col(v2));

            int rank1 = (index1.first - index1.second + 4) % 4;  // v1 -> v0
            int rank2 = (index2.first - index2.second + 4) % 4;  // v2 -> v0
            int orients[3] = {0};                                // == {0, 0, 0}
            if (v1 < v0) {
                variable_id[0] = -rshift90(variable_id[0], rank1);
                orients[0] = (rank1 + 2) % 4;
            } else {
                orients[0] = 0;
            }
            if (v2 < v1) {
                variable_id[1] = -rshift90(variable_id[1], rank2);
                orients[1] = (rank2 + 2) % 4;
            } else {
                variable_id[1] = rshift90(variable_id[1], rank1);
                orients[1] = rank1;
            }
            if (v2 < v0) {
                variable_id[2] = rshift90(variable_id[2], rank2);
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
        entities::DisjointOrientTree disajoint_orient_tree = entities::DisjointOrientTree(F.cols());
        // merge the whole face graph except for the singularity in which there exists a spanning tree
        // which contains consistent orientation
        std::vector<int> sharpUE(E2D.size());
        for (int i = 0; i < m_sharp_edges.size(); ++i) {
            if (m_sharp_edges[i]) {
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
        }
        {
            std::vector<int> total_flows(num_sharp_component);
            // check if each component is full-flow
            for (int i = 0; i < m_face_edge_ids.size(); ++i) {
                Vector2i diff(0, 0);
                for (int j = 0; j < 3; ++j) {
                    int orient = m_face_edge_orientation[i][j];
                    diff += rshift90(m_edge_difference[m_face_edge_ids[i][j]], orient);
                }
                total_flows[sharp_colors[i]] += diff[0] + diff[1];
            }

            // build "variable"
            m_variables.resize(m_edge_difference.size() * 2, std::make_pair(Vector2i(-1, -1), 0));
            for (int i = 0; i < m_face_edge_ids.size(); ++i) {
                for (int j = 0; j < 3; ++j) {
                    Vector2i sign = rshift90(Vector2i(1, 1), m_face_edge_orientation[i][j]);
                    int eid = m_face_edge_ids[i][j];
                    Vector2i index = rshift90(Vector2i(eid * 2, eid * 2 + 1), m_face_edge_orientation[i][j]);
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
            std::vector<std::vector<std::pair<int, int>>> modified_variables[2];
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
                Vector2i index = rshift90(Vector2i(e * 2 + 1, e * 2 + 2), m_face_edge_orientation[i][j]);
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
        std::vector<std::pair<Vector2i, int>> arcs;
        std::vector<int> arc_ids;
        entities::DisjointTree tree(m_face_edge_ids.size() * 2);
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
                diff += rshift90(m_edge_difference[m_face_edge_ids[i][j]], orient);
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
                Vector2i sign = rshift90(Vector2i(1, 1), m_face_edge_orientation[i][j]);
                int eid = m_face_edge_ids[i][j];
                Vector2i index = rshift90(Vector2i(eid * 2, eid * 2 + 1), m_face_edge_orientation[i][j]);
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
        std::vector<std::vector<std::pair<int, int>>> modified_variables[2];
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
            std::vector<entities::DEdge> &edge_values,
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
                diffs[eid] = rshift90(edge_diff[face_edgeIds[i][j]], face_edgeOrients[i][j]);
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
                auto value = compat_orientation_extrinsic_index_4(
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
                auto value = compat_orientation_extrinsic_index_4(
                        Q.col(v), N.col(v), face_spaces[f0].q, face_spaces[f0].n);
                if (F(j, f0) != v) orient += 2;
                face_edgeOrients[f0][j] = (orient + value.second - value.first + 4) % 4;
            }
            face_spaces[f0].d = d;
            for (int j = 0; j < 3; ++j) {
                int eid = face_edgeIds[f0][j];
                int orient = face_edgeOrients[f0][j];
                auto diff = rshift90(diffs[f0 * 3 + j], (4 - orient) % 4);
                edge_diff[eid] = diff;
            }
        };
        auto FixOrient = [&](int f0) {
            for (int j = 0; j < 3; ++j) {
                auto diff = edge_diff[face_edgeIds[f0][j]];
                if (rshift90(diff, face_edgeOrients[f0][j]) != diffs[f0 * 3 + j]) {
                    int orient = 0;
                    while (orient < 4 && rshift90(diff, orient) != diffs[f0 * 3 + j]) orient += 1;
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
            edge_values[eid0] = entities::DEdge(v0, vn);

            eid1 = edge_values.size();
            sharp_eid1 = sharp_eid;
            edge_values.push_back(entities::DEdge(vn, v1));
            edge_diff.push_back(Vector2i());

            eid0p = edge_values.size();
            sharp_eid0p = 0;
            edge_values.push_back(entities::DEdge(vn, v0p));
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
                while (rshift90(D01, orient) != Ds10) orient += 1;
                Vector2i Dsn0 = rshift90(D0n, orient);

                F.col(f1) << vn, v0, v1p;
                eid1p = edge_values.size();
                sharp_eid1p = 0;
                edge_values.push_back(entities::DEdge(vn, v1p));
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

            const int e0p = E2E[dedge_prev_3(e0)], e0n = E2E[dedge_next_3(e0)];

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
                const int e1p = E2E[dedge_prev_3(e1)], e1n = E2E[dedge_next_3(e1)];
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

    void Parametrizer::compute_index_map(Hierarchy &hierarchy, int with_scale) {
        auto &V = hierarchy.m_vertices[0];
        auto &F = hierarchy.m_faces;
        auto &Q = hierarchy.m_orientation[0];
        auto &N = hierarchy.m_normals[0];
        auto &O = hierarchy.m_positions[0];
        auto &S = hierarchy.m_scales[0];

        build_edge_info();

        // Constraints for the integer optimization

        for (int i = 0; i < m_sharp_edges.size(); ++i) {
            if (m_sharp_edges[i]) {
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
        std::map<int, std::pair<Vector3d, Vector3d>> sharp_constraints;
        std::set<int> sharpvert;
        for (int i = 0; i < m_sharp_edges.size(); ++i) {
            if (m_sharp_edges[i]) {
                sharpvert.insert(F(i % 3, i / 3));
                sharpvert.insert(F((i + 1) % 3, i / 3));
            }
        }

        m_allow_changes.resize(m_edge_difference.size() * 2, 1);
        for (int i = 0; i < m_sharp_edges.size(); ++i) {
            int e = m_face_edge_ids[i / 3][i % 3];
            if (sharpvert.count(m_edge_values[e].x) && sharpvert.count(m_edge_values[e].y)) {
                if (m_sharp_edges[i] != 0) {
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
                m_face_edge_ids, m_sharp_edges,
                m_singularities, 1
        );

        m_allow_changes.clear();
        m_allow_changes.resize(m_edge_difference.size() * 2, 1);
        for (int i = 0; i < m_sharp_edges.size(); ++i) {
            if (m_sharp_edges[i] == 0) continue;
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
                m_face_edge_ids, m_sharp_edges,
                m_singularities, 1
        );

        std::set<int> sharp_vertices;
        for (int i = 0; i < m_sharp_edges.size(); ++i) {
            if (m_sharp_edges[i] == 1) {
                sharp_vertices.insert(F(i % 3, i / 3));
                sharp_vertices.insert(F((i + 1) % 3, i / 3));
            }
        }

        Optimizer::optimize_positions_sharp(
                hierarchy,
                m_edge_values,
                m_edge_difference,
                m_sharp_edges,
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
        std::map<int, std::pair<Vector3d, Vector3d>> compact_sharp_constraints;
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

        std::vector<std::vector<int>> v2o(V.cols());
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
                            auto index = compat_orientation_extrinsic_index_4(q_1, n_1, q_2, n_2);
                            double s_x1 = S(0, v1), s_y1 = S(1, v1);
                            double s_x2 = S(0, v2), s_y2 = S(1, v2);
                            int rank_diff = (index.second + 4 - index.first) % 4;
                            if (rank_diff % 2 == 1) std::swap(s_x2, s_y2);
                            Vector3d qd_x = 0.5 * (rotate90_by(q_2, n_2, rank_diff) + q_1);
                            Vector3d qd_y = 0.5 * (rotate90_by(q_2_y, n_2, rank_diff) + q_1_y);
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
            Vector2i d1 = rshift90(m_edge_difference[m_face_edge_ids[i][0]], m_face_edge_orientation[i][0]);
            Vector2i d2 = rshift90(m_edge_difference[m_face_edge_ids[i][1]], m_face_edge_orientation[i][1]);
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
                compute_direct_graph_quad(m_positions_compact, m_faces_compact, m_V2E_compact, m_E2E_compact,
                                          m_boundary_compact, m_non_manifold_compact);
            } else {
                break;
            }
        }
        std::vector<std::vector<int>> v_dedges(m_V2E_compact.size());
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
        compute_direct_graph_quad(m_positions_compact, m_faces_compact, m_V2E_compact, m_E2E_compact,
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
            std::priority_queue<std::pair<int, int>> prior_queue;
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
                compute_direct_graph_quad(m_positions_compact, m_faces_compact, m_V2E_compact, m_E2E_compact,
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
        compute_direct_graph_quad(m_positions_compact, m_faces_compact, m_V2E_compact, m_E2E_compact,
                                  m_boundary_compact,
                                  m_non_manifold_compact);
        {
            compute_direct_graph_quad(m_positions_compact, m_faces_compact, m_V2E_compact, m_E2E_compact,
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
            std::vector<std::vector<int>> v_dedges(m_V2E_compact.size());
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

    void Parametrizer::close_hole(std::vector<int> &loop_vertices) {
        std::vector<std::vector<int>> loop_vertices_array;
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
                std::vector<int> vertices = {loop_vertices[0], loop_vertices[1], loop_vertices[seg1],
                                             loop_vertices[seg2]};
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
        disajoint_tree = entities::DisjointTree(V.cols());
        auto &diffs = fh.mEdgeDiff.front();
        for (int i = 0; i < diffs.size(); ++i) {
            if (diffs[i] == Vector2i::Zero()) {
                disajoint_tree.Merge(m_edge_values[i].x, m_edge_values[i].y);
            }
        }
        disajoint_tree.BuildCompactParent();
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
        int num_v = disajoint_tree.CompactNum();
        m_vertices_set.resize(num_v);
        m_positions_compact.resize(num_v, Vector3d::Zero());
        m_orientations_compact.resize(num_v, Vector3d::Zero());
        m_normals_compact.resize(num_v, Vector3d::Zero());
        m_counter.resize(num_v, 0);
        for (int i = 0; i < O.cols(); ++i) {
            int compact_v = disajoint_tree.Index(i);
            m_vertices_set[compact_v].push_back(i);
            m_positions_compact[compact_v] += O.col(i);
            m_normals_compact[compact_v] = m_normals_compact[compact_v] * m_counter[compact_v] + N.col(i);
            m_normals_compact[compact_v].normalize();
            if (m_counter[compact_v] == 0)
                m_orientations_compact[compact_v] = Q.col(i);
            else {
                auto pairs = compat_orientation_extrinsic_4(m_orientations_compact[compact_v],
                                                            m_normals_compact[compact_v],
                                                            Q.col(i), N.col(i));
                m_orientations_compact[compact_v] = (pairs.first * m_counter[compact_v] + pairs.second).normalized();
            }
            m_counter[compact_v] += 1;
        }
        for (int i = 0; i < m_positions_compact.size(); ++i) {
            m_positions_compact[i] /= m_counter[i];
        }

        build_triangle_manifold(disajoint_tree, edge, face, m_edge_values, F2E, E2F, EdgeDiff, FQ);
    }

    void Parametrizer::build_triangle_manifold(
            entities::DisjointTree &disajoint_tree,
            std::vector<int> &edge,
            std::vector<int> &face,
            std::vector<entities::DEdge> &edge_values,
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
        std::vector<std::vector<int>> Vs;
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
            std::vector<std::vector<int>> vert_to_dedge(num_v);
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
            Vector2i d1 = rshift90(EdgeDiff[triangle_edges[i][0]], triangle_orients[i][0]);
            Vector2i d2 = rshift90(EdgeDiff[triangle_edges[i][1]], triangle_orients[i][1]);
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
        compute_direct_graph(NV, NF, NV2E, NE2E, NB, NN);

        std::map<entities::DEdge, std::pair<Vector3i, Vector3i>> quads;
        for (int i = 0; i < triangle_vertices.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                int e = triangle_edges[i][j];
                int v1 = triangle_vertices[i][j];
                int v2 = triangle_vertices[i][(j + 1) % 3];
                int v3 = triangle_vertices[i][(j + 2) % 3];
                if (abs(EdgeDiff[e][0]) == 1 && abs(EdgeDiff[e][1]) == 1) {
                    entities::DEdge edge(v1, v2);
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

        compute_direct_graph_quad(
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
            compute_direct_graph_quad(
                    m_positions_compact,
                    m_faces_compact,
                    m_V2E_compact,
                    m_E2E_compact,
                    m_boundary_compact,
                    m_non_manifold_compact
            );
        }
        find_fix_holes();
        compute_direct_graph_quad(
                m_positions_compact,
                m_faces_compact,
                m_V2E_compact,
                m_E2E_compact,
                m_boundary_compact,
                m_non_manifold_compact
        );
    }

    // scale
    void Parametrizer::estimate_slope() {
        auto &mF = m_hierarchy.m_faces;
        auto &mQ = m_hierarchy.m_orientation[0];
        auto &mN = m_hierarchy.m_normals[0];
        auto &mV = m_hierarchy.m_vertices[0];
        m_faces_slope.resize(2, mF.cols());
        m_faces_orientation.resize(3, mF.cols());
        for (int i = 0; i < mF.cols(); ++i) {
            const Vector3d &n = m_faces_normals.col(i);
            const Vector3d &q_1 = mQ.col(mF(0, i)), &q_2 = mQ.col(mF(1, i)), &q_3 = mQ.col(mF(2, i));
            const Vector3d &n_1 = mN.col(mF(0, i)), &n_2 = mN.col(mF(1, i)), &n_3 = mN.col(mF(2, i));
            Vector3d q_1n = rotate_vector_into_plane(q_1, n_1, n);
            Vector3d q_2n = rotate_vector_into_plane(q_2, n_2, n);
            Vector3d q_3n = rotate_vector_into_plane(q_3, n_3, n);

            auto p = compat_orientation_extrinsic_4(q_1n, n, q_2n, n);
            Vector3d q = (p.first + p.second).normalized();
            p = compat_orientation_extrinsic_4(q, n, q_3n, n);
            q = (p.first * 2 + p.second);
            q = q - n * q.dot(n);
            m_faces_orientation.col(i) = q.normalized();
        }
        for (int i = 0; i < mF.cols(); ++i) {
            double step = m_hierarchy.m_scale * 1.f;

            const Vector3d &n = m_faces_normals.col(i);
            Vector3d p = (mV.col(mF(0, i)) + mV.col(mF(1, i)) + mV.col(mF(2, i))) * (1.0 / 3.0);
            Vector3d q_x = m_faces_orientation.col(i), q_y = n.cross(q_x);
            Vector3d q_xl = -q_x, q_xr = q_x;
            Vector3d q_yl = -q_y, q_yr = q_y;
            Vector3d q_yl_unfold = q_y, q_yr_unfold = q_y, q_xl_unfold = q_x, q_xr_unfold = q_x;
            int f;
            double tx, ty, len;

            f = i;
            len = step;
            TravelField(p, q_xl, len, f, m_hierarchy.m_E2E, mV, mF, m_faces_normals, m_faces_orientation, mQ, mN,
                        m_triangle_space, &tx, &ty,
                        &q_yl_unfold);

            f = i;
            len = step;
            TravelField(p, q_xr, len, f, m_hierarchy.m_E2E, mV, mF, m_faces_normals, m_faces_orientation, mQ, mN,
                        m_triangle_space, &tx, &ty,
                        &q_yr_unfold);

            f = i;
            len = step;
            TravelField(p, q_yl, len, f, m_hierarchy.m_E2E, mV, mF, m_faces_normals, m_faces_orientation, mQ, mN,
                        m_triangle_space, &tx, &ty,
                        &q_xl_unfold);

            f = i;
            len = step;
            TravelField(p, q_yr, len, f, m_hierarchy.m_E2E, mV, mF, m_faces_normals, m_faces_orientation, mQ, mN,
                        m_triangle_space, &tx, &ty,
                        &q_xr_unfold);
            double dSx = (q_yr_unfold - q_yl_unfold).dot(q_x) / (2.0f * step);
            double dSy = (q_xr_unfold - q_xl_unfold).dot(q_y) / (2.0f * step);
            m_faces_slope.col(i) = Vector2d(dSx, dSy);
        }

        std::vector<double> areas(mV.cols(), 0.0);
        for (int i = 0; i < mF.cols(); ++i) {
            Vector3d p1 = mV.col(mF(1, i)) - mV.col(mF(0, i));
            Vector3d p2 = mV.col(mF(2, i)) - mV.col(mF(0, i));
            double area = p1.cross(p2).norm();
            for (int j = 0; j < 3; ++j) {
                auto index = compat_orientation_extrinsic_index_4(m_faces_orientation.col(i), m_faces_normals.col(i),
                                                                  mQ.col(mF(j, i)),
                                                                  mN.col(mF(j, i)));
                double scaleX = m_faces_slope.col(i).x(), scaleY = m_faces_slope.col(i).y();
                if (index.first != index.second % 2) {
                    std::swap(scaleX, scaleY);
                }
                if (index.second >= 2) {
                    scaleX = -scaleX;
                    scaleY = -scaleY;
                }
                m_hierarchy.m_areas[0].col(mF(j, i)) += area * Vector2d(scaleX, scaleY);
                areas[mF(j, i)] += area;
            }
        }
        for (int i = 0; i < mV.cols(); ++i) {
            if (areas[i] != 0)
                m_hierarchy.m_areas[0].col(i) /= areas[i];
        }
        for (int l = 0; l < m_hierarchy.m_areas.size() - 1; ++l) {
            const MatrixXd &K = m_hierarchy.m_areas[l];
            MatrixXd &K_next = m_hierarchy.m_areas[l + 1];
            auto &toUpper = m_hierarchy.mToUpper[l];
            for (int i = 0; i < toUpper.cols(); ++i) {
                Vector2i upper = toUpper.col(i);
                Vector2d k0 = K.col(upper[0]);

                if (upper[1] != -1) {
                    Vector2d k1 = K.col(upper[1]);
                    k0 = 0.5 * (k0 + k1);
                }

                K_next.col(i) = k0;
            }
        }
    }
}
