#include <atomic>
#include <fstream>
#include <vector>

#include "spdlog/spdlog.h"
#include "mathext.h"

namespace mathext {
    class SVDSolver {
        std::tuple<float, float> rotate_xy(
            float x,
            float y,
            float cosine,
            float sine
        ) {
            const float x_new = cosine * x - sine * y;
            const float y_new = sine * x + cosine * y;
            return std::make_tuple(x_new, y_new);
        }

        std::tuple<float, float> rotate_q_xy(
            float x,
            float y,
            float a,
            float cosine,
            float sine
        ) {
            const float cc = cosine * cosine;
            const float ss = sine * sine;
            const float mx = 2.0f * cosine * sine * a;
            const float x_new = cc * x - mx + ss * y;
            const float y_new = ss * x + mx + cc * y;
            return std::make_tuple(x_new, y_new);
        }

        std::tuple<float, float> givens_rotation_coefficients(
            const float pp,
            const float pq,
            const float qq
        ) {
            const float tau = (qq - pp) / (2.0f * pq + 1e-6f);
            const float stt = std::sqrtf(1.0 + tau * tau);
            const float t = 1.0f / ((tau >= 0.0f) ? (tau + stt) : (tau - stt));
            const float cosine = 1.0f / std::sqrt(1.0f + t * t);
            const float sine = t * cosine;
            return std::make_tuple(cosine, sine);
        }

        void rotate01(
            Matrix3f &vtav,
            Matrix3f &v
        ) {
            if (vtav(0, 1) == 0.0f) return;

            const auto [cosine, sine] = givens_rotation_coefficients(vtav(0, 0), vtav(0, 1), vtav(1, 1));
            const auto [qx, qy] = rotate_q_xy(vtav(0, 0), vtav(1, 1), vtav(0, 1), cosine, sine);
            const auto [x, y] = rotate_xy(vtav(0, 2), vtav(1, 2), cosine, sine);

            vtav(0, 0) = qx;
            vtav(1, 1) = qy;
            vtav(0, 2) = x;
            vtav(1, 2) = y;
            vtav(0, 1) = 0.0f;

            for (int i = 0; i < 3; ++i) {
                const auto [x, y] = rotate_xy(v(i, 0), v(i, 1), cosine, sine);
                v(i, 0) = x;
                v(i, 1) = y;
            }
        }

        void rotate02(
            Matrix3f &vtav,
            Matrix3f &v
        ) {
            if (vtav(0, 2) == 0.0f) return;

            const auto [cosine, sine] = givens_rotation_coefficients(vtav(0, 0), vtav(0, 2), vtav(2, 2));
            const auto [qx, qy] = rotate_q_xy(vtav(0, 0), vtav(2, 2), vtav(0, 2), cosine, sine);
            const auto [x, y] = rotate_xy(vtav(0, 1), vtav(1, 2), cosine, sine);

            vtav(0, 0) = qx;
            vtav(2, 2) = qy;
            vtav(0, 1) = x;
            vtav(1, 2) = y;
            vtav(0, 2) = 0.0f;

            for (int i = 0; i < 3; ++i) {
                const auto [x, y] = rotate_xy(v(i, 0), v(i, 2), cosine, sine);
                v(i, 0) = x;
                v(i, 2) = y;
            }
        }

        void rotate12(
            Matrix3f &vtav,
            Matrix3f &v
        ) {
            if (vtav(1, 2) == 0.0f) return;

            const auto [cosine, sine] = givens_rotation_coefficients(vtav(1, 1), vtav(1, 2), vtav(2, 2));
            const auto [qx, qy] = rotate_q_xy(vtav(1, 1), vtav(2, 2), vtav(1, 2), cosine, sine);
            const auto [x, y] = rotate_xy(vtav(0, 1), vtav(0, 2), cosine, sine);

            vtav(1, 1) = qx;
            vtav(2, 2) = qy;
            vtav(0, 1) = x;
            vtav(0, 2) = y;
            vtav(1, 2) = 0.0f;

            for (int i = 0; i < 3; ++i) {
                const auto [x, y] = rotate_xy(v(i, 1), v(i, 2), cosine, sine);
                v(i, 1) = x;
                v(i, 2) = y;
            }
        }


        std::tuple<Vector3f, Matrix3f> solve_symmetric(
            const Matrix3f &ATA,
            const int iterations = 10
        ) {
            Matrix3f vtav = ATA;
            Matrix3f v = MatrixXf::Identity(3, 3);
            for (int i = 0; i < iterations; ++i) {
                rotate01(vtav, v);
                rotate02(vtav, v);
                rotate12(vtav, v);
            }
            const Vector3f sigma = {vtav(0, 0), vtav(1, 1), vtav(2, 2)};
            return std::make_tuple(sigma, v);
        }

        float inverse_determinant(
            const float x,
            const float tolerance = 1e-1f
        ) {
            if (std::abs(x) < tolerance || std::abs(1.0f / x) < tolerance) return 0.0f;
            return 1.0f / x;
        }

        Matrix3f pseudoinverse(
            const Vector3f &sigma,
            const Matrix3f &V
        ) {
            const auto determinant = Vector3f(
                inverse_determinant(sigma[0]),
                inverse_determinant(sigma[1]),
                inverse_determinant(sigma[2])
            );

            const Matrix3f a = V.array().rowwise() * determinant.array().transpose();
            return a * V.transpose();
        }

    public:
        /**
         * Solves the linear system ATA * x = ATb using SVD.
         * @param ATA The matrix A^T * A.
         * @param ATb The vector A^T * b.
         * @return The solution vector x.
         */
        Vector3f solve(
            const Matrix3f &ATA,
            const Vector3f &ATb
        ) {
            const auto [sigma, V] = solve_symmetric(ATA);
            const auto V_inv = pseudoinverse(sigma, V);

            return V_inv * ATb;
        }
    };


    MatrixXf clip(
        const MatrixXf &mat,
        float minVal,
        float maxVal
    ) {
        return mat.unaryExpr([minVal, maxVal](float val) {
            return std::min(std::max(val, minVal), maxVal);
        });
    }

    float frobenius_norm_off_diagonal(const Eigen::MatrixXf &A) {
        int n = A.rows();
        float sum = 0.0f;

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += A(i, j) * A(i, j);
                }
            }
        }

        return std::sqrt(sum);
    }

    float percentile(const Eigen::VectorXf &vec, float percentile) {
        // Copy the data from the Eigen vector to a standard vector for sorting
        Eigen::VectorXf vec_eval = vec.eval();
        std::vector<float> sorted_vec(vec_eval.data(), vec_eval.data() + vec_eval.size());

        // Sort the data
        std::sort(sorted_vec.begin(), sorted_vec.end());

        // Compute the percentile index
        float percentile_index = percentile * (sorted_vec.size() - 1);
        int lower_index = static_cast<int>(percentile_index);
        int upper_index = lower_index + 1;

        // Get the lower and upper values
        float lower_value = sorted_vec[lower_index];
        float upper_value = (upper_index < sorted_vec.size()) ? sorted_vec[upper_index] : lower_value;

        // Linear interpolation if needed
        float weight = percentile_index - lower_index;
        float percentile_value = lower_value * (1.0f - weight) + upper_value * weight;

        return percentile_value;
    }

    // Mesh math
    Vector3f face_centroid(
        entities::Mesh &mesh,
        const entities::Mesh::FaceHandle &face
    ) {
        int total = 0;
        Vector3f centroid = Vector3f::Zero();
        for (auto it_face_vertex = mesh.fv_iter(face);
             it_face_vertex.is_valid(); ++it_face_vertex) {
            const auto p = mesh.point(*it_face_vertex);
            centroid += Vector3f(p[0], p[1], p[2]);
            total++;
        }

        centroid /= total;

        return centroid;
    }

    MatrixXf face_centroids_ring(
        entities::Mesh &mesh,
        const entities::Mesh::VertexHandle vertex
    ) {
        std::vector<Vector3f> centroids_list;
        for (auto it_face = mesh.vf_iter(vertex);
             it_face.is_valid(); ++it_face) {
            const auto centroid = face_centroid(mesh, *it_face);
            centroids_list.push_back(centroid);
        }

        MatrixXf centroids(centroids_list.size(), 3);
        for (int i = 0; i < centroids_list.size(); ++i) {
            centroids.row(i) = centroids_list[i];
        }

        return centroids;
    }

    MatrixXf count_unique(
        const MatrixXf &mat
    ) {
        // Create an unordered map to store unique elements and their counts
        std::map<float, int> elementCounts;

        // Iterate over all elements of the matrix and populate the map with counts
        for (int i = 0; i < mat.rows(); ++i) {
            for (int j = 0; j < mat.cols(); ++j) {
                elementCounts[mat(i, j)]++;
            }
        }

        // Create a MatrixXf with two columns: one for unique elements, one for their counts
        int uniqueCount = elementCounts.size();
        MatrixXf uniqueElementsWithCounts(uniqueCount, 2);

        // Fill the matrix with unique elements and their counts
        int index = 0;
        for (const auto &[element, count]: elementCounts) {
            uniqueElementsWithCounts(index, 0) = element; // Unique element
            uniqueElementsWithCounts(index, 1) = count; // Count of that element
            index++;
        }

        return uniqueElementsWithCounts;
    }

    Vector4f _intersect_planes(
        const MatrixXf &vertices,
        const MatrixXf &normals
    ) {
        const auto svd_vmul_sym = [](Matrix3f a, Vector3f v) {
            return Vector3f{
                a.row(0).dot(v),
                a(0, 1) * v[0] + a(1, 1) * v[1] + a(1, 2) * v[2],
                a(0, 2) * v[0] + a(1, 2) * v[1] + a(2, 2) * v[2]
            };
        };

        const Matrix3f ATA = (normals.transpose() * normals).triangularView<Upper>();
        const VectorXf b = (vertices.array() * normals.array()).rowwise().sum();
        Vector3f ATb = (b.asDiagonal() * normals).colwise().sum();

        const Vector3f center = vertices.colwise().mean();
        const Vector3f constraint_center = svd_vmul_sym(ATA, center);
        ATb -= constraint_center.transpose();

        Vector3f x = SVDSolver().solve(ATA, ATb);

        const Vector3f differences = ATb - svd_vmul_sym(ATA, x);
        const float error = differences.dot(differences);

        x += center;

        return {x[0], x[1], x[2], error};
    }

    Vector4f intersect_planes(
        const MatrixXf &vertices,
        const MatrixXf &normals,
        const int patience
    ) {
        assert(vertices.rows() == normals.rows());

        const float max_displacement = (vertices.rowwise() - vertices.colwise().mean())
                .rowwise()
                .norm()
                .minCoeff();

        float displacement = 0;
        Vector4f solution;
        int iteration = 0;
        do {
            solution = _intersect_planes(vertices, normals);
            displacement = (solution.head<3>()).norm();

            iteration++;
        } while (displacement > max_displacement && iteration < patience);

        return solution;
    }


    std::tuple<MatrixXf, float, Vector3f> normalize(
        const MatrixXf &vertices
    ) {
        const Vector3f upper = vertices.colwise().maxCoeff();
        const Vector3f lower = vertices.colwise().minCoeff();

        const float scale = (upper - lower).maxCoeff() * 0.5;
        const Vector3f offset = (upper + lower) * 0.5;
        const MatrixXf vertices_scaled = (vertices.rowwise() - offset.transpose()) / scale;

        return std::make_tuple(
            vertices_scaled,
            scale,
            offset
        );
    }

    MatrixXf denormalize(
        const MatrixXf &vertices,
        const float scale,
        const Vector3f &offset
    ) {
        return (vertices * scale).rowwise() + offset.transpose();
    }

    std::vector<std::vector<entities::Mesh::VertexHandle>> find_boundary_vertices(const entities::Mesh &mesh) {
        std::vector<std::vector<entities::Mesh::VertexHandle> > boundaries;

        for (auto heh: mesh.halfedges()) {
            if (mesh.is_boundary(heh)) {
                OpenMesh::HalfedgeHandle start = heh;
                std::vector<entities::Mesh::VertexHandle> boundary;
                do {
                    boundary.emplace_back(mesh.to_vertex_handle(heh));
                    heh = mesh.next_halfedge_handle(heh);
                } while (heh != start);
                boundaries.emplace_back(boundary);
            }
        }

        return boundaries;
    }
}

namespace mathext {
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
        std::vector<std::pair<uint32_t, uint32_t> > tmp(F.size());

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
        std::vector<std::pair<uint32_t, uint32_t> > tmp(F.size() * deg);

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
