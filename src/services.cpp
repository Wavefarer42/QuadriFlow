#include <fstream>

#include "spdlog/spdlog.h"

#include "services.h"
#include "field-math.h"

namespace services {

    entities::QuadMesh MeshService::load_trimesh_from_file(const std::string &filename) const {
        spdlog::info("Loading triangle mesh from file from {}", filename);

        const auto mesh = this->mesh_dao.load_mesh_from_file(filename);

        // Validate that the mesh is a triangle mesh
        int n_triangles = 0;
        int n_non_triangles = 0;
        for (auto it_f = mesh.faces_begin(); it_f != mesh.faces_end(); ++it_f) {
            if (mesh.valence(*it_f) == 3) {
                n_triangles++;
            } else {
                n_non_triangles++;
            }
        }

        if (n_non_triangles > 0) {
            throw std::runtime_error(
                    std::format("Please provide a triangle mesh as input. Triangles={}, Non Triangle={}",
                                n_triangles, n_non_triangles)
            );
        }

        return mesh;
    }

    void MeshService::save_quadmesh_to_file(const std::string &filename, Parametrizer &field) const {
        spdlog::info("Saving mesh to file {}", filename);

        std::ofstream os(filename);
        for (int i = 0; i < field.m_positions_compact.size(); ++i) {
            auto t = field.m_positions_compact[i] * field.m_normalize_scale + field.m_normalize_offset;
            os << "v " << t[0] << " " << t[1] << " " << t[2] << "\n";
        }
        for (int i = 0; i < field.m_faces_compact.size(); ++i) {
            os << "f " << field.m_faces_compact[i][0] + 1 << " " << field.m_faces_compact[i][1] + 1 << " "
               << field.m_faces_compact[i][2] + 1 << " " << field.m_faces_compact[i][3] + 1 << "\n";
        }
        os.close();
    }

    void MeshService::set_boundary_constraints(Hierarchy &hierarchy) const {
        spdlog::info("Setting boundary constraints");

        for (uint32_t i = 0; i < 3 * hierarchy.m_faces.cols(); ++i) {
            if (hierarchy.m_E2E[i] == -1) {
                uint32_t i0 = hierarchy.m_faces(i % 3, i / 3);
                uint32_t i1 = hierarchy.m_faces((i + 1) % 3, i / 3);
                Vector3d p0 = hierarchy.m_vertices[0].col(i0), p1 = hierarchy.m_vertices[0].col(i1);
                Vector3d edge = p1 - p0;
                if (edge.squaredNorm() > 0) {
                    edge.normalize();
                    hierarchy.m_position_constraints[0].col(i0) = p0;
                    hierarchy.m_position_constraints[0].col(i1) = p1;
                    hierarchy.m_orientation_constraint[0].col(i0) = hierarchy.m_orientation_constraint[0].col(
                            i1) = edge;
                    hierarchy.m_orientation_constraint_weight[0][i0] = hierarchy.m_orientation_constraint_weight[0][i1] = hierarchy.m_position_constraint_weights[0][i0] = hierarchy.m_position_constraint_weights[0][i1] =
                            1.0;
                }
            }
        }

        hierarchy.propagateConstraints();
    }

    std::map<int, int> MeshService::find_orientation_singularities(Hierarchy &hierarchy) const {
        spdlog::info("Finding orientation singularities");

        const MatrixXd &normals = hierarchy.m_normals[0];
        const MatrixXi &faces = hierarchy.m_faces;
        MatrixXd &Q = hierarchy.m_orientation[0];

        std::map<int, int> singularities;
        for (int f = 0; f < faces.cols(); ++f) {
            int index = 0;
            int abs_index = 0;
            for (int k = 0; k < 3; ++k) {
                int i = faces(k, f), j = faces(k == 2 ? 0 : (k + 1), f);
                auto value = compat_orientation_extrinsic_index_4(
                        Q.col(i),
                        normals.col(i),
                        Q.col(j),
                        normals.col(j)
                );
                index += value.second - value.first;
                abs_index += std::abs(value.second - value.first);
            }
            int index_mod = modulo(index, 4);
            if (index_mod == 1 || index_mod == 3) {
                if (index >= 4 || index < 0) {
                    // TODO is the negative sign a marking?
                    Q.col(faces(0, f)) = -Q.col(faces(0, f));
                }
                singularities[f] = index_mod;
            }
        }

        return singularities;
    }

    std::tuple<std::map<int, Vector2i>, MatrixXi, MatrixXi> MeshService::find_position_singularities(
            Hierarchy &m_hierarchy,
            bool with_scale
    ) const {
        const MatrixXd &V = m_hierarchy.m_vertices[0];
        const MatrixXd &N = m_hierarchy.m_normals[0];
        const MatrixXd &Q = m_hierarchy.m_orientation[0];
        const MatrixXd &O = m_hierarchy.m_positions[0];
        const MatrixXi &F = m_hierarchy.m_faces;

        std::map<int, Vector2i> singularity_position;
        MatrixXi singularity_rank(F.rows(), F.cols());
        MatrixXi singularity_index(6, F.cols());

        for (int f = 0; f < F.cols(); ++f) {
            Vector2i index = Vector2i::Zero();
            uint32_t i0 = F(0, f), i1 = F(1, f), i2 = F(2, f);

            Vector3d q[3] = {Q.col(i0).normalized(), Q.col(i1).normalized(), Q.col(i2).normalized()};
            Vector3d n[3] = {N.col(i0), N.col(i1), N.col(i2)};
            Vector3d o[3] = {O.col(i0), O.col(i1), O.col(i2)};
            Vector3d v[3] = {V.col(i0), V.col(i1), V.col(i2)};

            int best[3];
            double best_dp = -std::numeric_limits<double>::infinity();
            for (int i = 0; i < 4; ++i) {
                Vector3d v0 = rotate90_by(q[0], n[0], i);
                for (int j = 0; j < 4; ++j) {
                    Vector3d v1 = rotate90_by(q[1], n[1], j);
                    for (int k = 0; k < 4; ++k) {
                        Vector3d v2 = rotate90_by(q[2], n[2], k);
                        double dp = std::min(std::min(v0.dot(v1), v1.dot(v2)), v2.dot(v0));
                        if (dp > best_dp) {
                            best_dp = dp;
                            best[0] = i;
                            best[1] = j;
                            best[2] = k;
                        }
                    }
                }
            }
            singularity_rank(0, f) = best[0];
            singularity_rank(1, f) = best[1];
            singularity_rank(2, f) = best[2];
            for (int k = 0; k < 3; ++k) q[k] = rotate90_by(q[k], n[k], best[k]);

            for (int k = 0; k < 3; ++k) {
                int kn = k == 2 ? 0 : (k + 1);
                double scale_x = m_hierarchy.m_scale, scale_y = m_hierarchy.m_scale,
                        scale_x_1 = m_hierarchy.m_scale, scale_y_1 = m_hierarchy.m_scale;
                if (with_scale) {
                    scale_x *= m_hierarchy.m_scales[0](0, F(k, f));
                    scale_y *= m_hierarchy.m_scales[0](1, F(k, f));
                    scale_x_1 *= m_hierarchy.m_scales[0](0, F(kn, f));
                    scale_y_1 *= m_hierarchy.m_scales[0](1, F(kn, f));
                    if (best[k] % 2 != 0) std::swap(scale_x, scale_y);
                    if (best[kn] % 2 != 0) std::swap(scale_x_1, scale_y_1);
                }
                double inv_scale_x = 1.0 / scale_x, inv_scale_y = 1.0 / scale_y,
                        inv_scale_x_1 = 1.0 / scale_x_1, inv_scale_y_1 = 1.0 / scale_y_1;
                std::pair<Vector2i, Vector2i> value = compat_position_extrinsic_index_4(
                        v[k], n[k], q[k], o[k], v[kn], n[kn], q[kn], o[kn], scale_x, scale_y, inv_scale_x,
                        inv_scale_y, scale_x_1, scale_y_1, inv_scale_x_1, inv_scale_y_1, nullptr);
                auto diff = value.first - value.second;
                index += diff;
                singularity_index(k * 2, f) = diff[0];
                singularity_index(k * 2 + 1, f) = diff[1];
            }

            if (index != Vector2i::Zero()) {
                singularity_position[f] = rshift90(index, best[0]);
            }
        }

        return std::make_tuple(singularity_position, singularity_rank, singularity_index);
    }

    std::tuple<MatrixXd, MatrixXd> MeshService::estimate_slope(
            Hierarchy &hierarchy,
            std::vector<MatrixXd> &triangle_space,
            MatrixXd &normals_faces
    ) const {
        spdlog::info("Estimating adaptive slope");

        auto &faces = hierarchy.m_faces;
        auto &orientation = hierarchy.m_orientation[0];
        auto &normals = hierarchy.m_normals[0];
        auto &vertices = hierarchy.m_vertices[0];

        MatrixXd faces_slope(2, faces.cols());
        MatrixXd faces_orientation(3, faces.cols());

        for (int i = 0; i < faces.cols(); ++i) {
            const Vector3d &n = normals_faces.col(i);
            const Vector3d &q_1 = orientation.col(faces(0, i)), &q_2 = orientation.col(
                    faces(1, i)), &q_3 = orientation.col(faces(2, i));
            const Vector3d &n_1 = normals.col(faces(0, i)), &n_2 = normals.col(faces(1, i)), &n_3 = normals.col(
                    faces(2, i));
            Vector3d q_1n = rotate_vector_into_plane(q_1, n_1, n);
            Vector3d q_2n = rotate_vector_into_plane(q_2, n_2, n);
            Vector3d q_3n = rotate_vector_into_plane(q_3, n_3, n);

            auto p = compat_orientation_extrinsic_4(q_1n, n, q_2n, n);
            Vector3d q = (p.first + p.second).normalized();
            p = compat_orientation_extrinsic_4(q, n, q_3n, n);
            q = (p.first * 2 + p.second);
            q = q - n * q.dot(n);
            faces_orientation.col(i) = q.normalized();
        }
        for (int i = 0; i < faces.cols(); ++i) {
            double step = hierarchy.m_scale * 1.f;

            const Vector3d &n = normals_faces.col(i);
            Vector3d p =
                    (vertices.col(faces(0, i)) + vertices.col(faces(1, i)) + vertices.col(faces(2, i))) * (1.0 / 3.0);
            Vector3d q_x = faces_orientation.col(i), q_y = n.cross(q_x);
            Vector3d q_xl = -q_x, q_xr = q_x;
            Vector3d q_yl = -q_y, q_yr = q_y;
            Vector3d q_yl_unfold = q_y, q_yr_unfold = q_y, q_xl_unfold = q_x, q_xr_unfold = q_x;
            int f;
            double tx, ty, len;

            f = i;
            len = step;
            TravelField(p, q_xl, len, f, hierarchy.m_E2E, vertices, faces, normals_faces, faces_orientation,
                        orientation, normals, triangle_space, &tx,
                        &ty,
                        &q_yl_unfold);

            f = i;
            len = step;
            TravelField(p, q_xr, len, f, hierarchy.m_E2E, vertices, faces, normals_faces, faces_orientation,
                        orientation, normals, triangle_space, &tx,
                        &ty,
                        &q_yr_unfold);

            f = i;
            len = step;
            TravelField(p, q_yl, len, f, hierarchy.m_E2E, vertices, faces, normals_faces, faces_orientation,
                        orientation, normals, triangle_space, &tx,
                        &ty,
                        &q_xl_unfold);

            f = i;
            len = step;
            TravelField(p, q_yr, len, f, hierarchy.m_E2E, vertices, faces, normals_faces, faces_orientation,
                        orientation, normals, triangle_space, &tx,
                        &ty,
                        &q_xr_unfold);
            double dSx = (q_yr_unfold - q_yl_unfold).dot(q_x) / (2.0f * step);
            double dSy = (q_xr_unfold - q_xl_unfold).dot(q_y) / (2.0f * step);
            faces_slope.col(i) = Vector2d(dSx, dSy);
        }

        std::vector<double> areas(vertices.cols(), 0.0);
        for (int i = 0; i < faces.cols(); ++i) {
            Vector3d p1 = vertices.col(faces(1, i)) - vertices.col(faces(0, i));
            Vector3d p2 = vertices.col(faces(2, i)) - vertices.col(faces(0, i));
            double area = p1.cross(p2).norm();
            for (int j = 0; j < 3; ++j) {
                auto index = compat_orientation_extrinsic_index_4(faces_orientation.col(i), normals_faces.col(i),
                                                                  orientation.col(faces(j, i)),
                                                                  normals.col(faces(j, i)));
                double scaleX = faces_slope.col(i).x(), scaleY = faces_slope.col(i).y();
                if (index.first != index.second % 2) {
                    std::swap(scaleX, scaleY);
                }
                if (index.second >= 2) {
                    scaleX = -scaleX;
                    scaleY = -scaleY;
                }
                hierarchy.m_areas[0].col(faces(j, i)) += area * Vector2d(scaleX, scaleY);
                areas[faces(j, i)] += area;
            }
        }
        for (int i = 0; i < vertices.cols(); ++i) {
            if (areas[i] != 0)
                hierarchy.m_areas[0].col(i) /= areas[i];
        }
        for (int l = 0; l < hierarchy.m_areas.size() - 1; ++l) {
            const MatrixXd &K = hierarchy.m_areas[l];
            MatrixXd &K_next = hierarchy.m_areas[l + 1];
            auto &toUpper = hierarchy.mToUpper[l];
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

        return std::make_tuple(faces_slope, faces_orientation);
    }

    entities::QuadMesh mesh_sdfn(const std::string path_model) {
        spdlog::info("Meshing SDFn via surface nets.");

        entities::QuadMesh mesh;
        return mesh;
    }

}