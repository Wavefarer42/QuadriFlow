#include "spdlog/spdlog.h"

#include "services.h"
#include "field-math.h"


namespace services {

    entities::QuadMesh MeshService::load_trimesh_from_file(const std::string &filename) {
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

    void MeshService::save_quadmesh_to_file(const std::string &filename, const entities::QuadMesh &mesh) {
        spdlog::info("Saving mesh to file {}", filename);

        this->mesh_dao.save_mesh_to_file(filename, mesh);
    }

    void MeshService::set_boundary_constraints(Hierarchy &hierarchy) {
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

    std::map<int, int> MeshService::find_orientation_singularities(Hierarchy &hierarchy) {
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

    std::tuple<MatrixXd, MatrixXd> MeshService::estimate_slope(
            Hierarchy &hierarchy,
            std::vector<MatrixXd> &triangle_space,
            MatrixXd &normals_faces
    ) {
        spdlog::info("Estimating adaptive slope");

        auto &faces = hierarchy.m_faces;
        auto &orientation = hierarchy.m_orientation[0];
        auto &normals = hierarchy.m_normals[0];
        auto &vertices = hierarchy.m_vertices[0];

        MatrixXd faces_slope(2, faces.cols());
        MatrixXd faces_orientation(3, faces.cols());

        for (int i = 0; i < faces.cols(); ++i) {
            const Vector3d &n = normals_faces.col(i);
            const Vector3d &q_1 = orientation.col(faces(0, i)), &q_2 = orientation.col(faces(1, i)), &q_3 = orientation.col(faces(2, i));
            const Vector3d &n_1 = normals.col(faces(0, i)), &n_2 = normals.col(faces(1, i)), &n_3 = normals.col(faces(2, i));
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
            Vector3d p = (vertices.col(faces(0, i)) + vertices.col(faces(1, i)) + vertices.col(faces(2, i))) * (1.0 / 3.0);
            Vector3d q_x = faces_orientation.col(i), q_y = n.cross(q_x);
            Vector3d q_xl = -q_x, q_xr = q_x;
            Vector3d q_yl = -q_y, q_yr = q_y;
            Vector3d q_yl_unfold = q_y, q_yr_unfold = q_y, q_xl_unfold = q_x, q_xr_unfold = q_x;
            int f;
            double tx, ty, len;

            f = i;
            len = step;
            TravelField(p, q_xl, len, f, hierarchy.m_E2E, vertices, faces, normals_faces, faces_orientation, orientation, normals, triangle_space, &tx,
                        &ty,
                        &q_yl_unfold);

            f = i;
            len = step;
            TravelField(p, q_xr, len, f, hierarchy.m_E2E, vertices, faces, normals_faces, faces_orientation, orientation, normals, triangle_space, &tx,
                        &ty,
                        &q_yr_unfold);

            f = i;
            len = step;
            TravelField(p, q_yl, len, f, hierarchy.m_E2E, vertices, faces, normals_faces, faces_orientation, orientation, normals, triangle_space, &tx,
                        &ty,
                        &q_xl_unfold);

            f = i;
            len = step;
            TravelField(p, q_yr, len, f, hierarchy.m_E2E, vertices, faces, normals_faces, faces_orientation, orientation, normals, triangle_space, &tx,
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
                auto index = compat_orientation_extrinsic_index_4(faces_orientation.col(i), normals_faces.col(i), orientation.col(faces(j, i)),
                                                                  normals.col(faces(j, i)));
                double scaleX = faces_slope.col(i).x(), scaleY = faces_slope.col(i).y();
                if (index.first != index.second % 2) {
                    std::swap(scaleX, scaleY);
                }
                if (index.second >= 2) {
                    scaleX = -scaleX;
                    scaleY = -scaleY;
                }
                hierarchy.mK[0].col(faces(j, i)) += area * Vector2d(scaleX, scaleY);
                areas[faces(j, i)] += area;
            }
        }
        for (int i = 0; i < vertices.cols(); ++i) {
            if (areas[i] != 0)
                hierarchy.mK[0].col(i) /= areas[i];
        }
        for (int l = 0; l < hierarchy.mK.size() - 1; ++l) {
            const MatrixXd &K = hierarchy.mK[l];
            MatrixXd &K_next = hierarchy.mK[l + 1];
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
}