#include <fstream>
#include <algorithm>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

#include "services.h"
#include "field-math.h"
#include "optimizer.h"
#include "adapters.h"
#include "surfacenets.h"
#include "transformations.h"
#include "sdfn.h"


template<typename Derived>
Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>
clip(const Eigen::MatrixBase<Derived> &mat, typename Derived::Scalar minVal, typename Derived::Scalar maxVal) {
    return mat.unaryExpr([minVal, maxVal](typename Derived::Scalar val) {
        return std::min(std::max(val, minVal), maxVal);
    });
}

template<typename MatrixType>
typename MatrixType::Scalar frobenius_norm_off_diagonal(const MatrixType &A) {
    using Scalar = typename MatrixType::Scalar;

    int n = A.rows();
    Scalar sum = 0;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                sum += A(i, j) * A(i, j);
            }
        }
    }

    return std::sqrt(sum);
}

// Mesh functions
Vector3f face_centroid(
        entities::Mesh &mesh,
        entities::Mesh::VertexFaceIter &face
) {

    int total = 0;
    Vector3f centroid = Vector3f::Zero();
    for (auto it_face_vertex = mesh.fv_iter(*face);
         it_face_vertex.is_valid(); ++it_face_vertex) {
        const auto p = mesh.point(*it_face_vertex);
        centroid += Vector3f(p[0], p[1], p[2]);
        total++;
    }

    centroid /= total;

    return centroid;
};

MatrixXf face_centroids_ring(
        entities::Mesh &mesh,
        entities::Mesh::VertexHandle vertex
) {
    std::cout << "Face: ";
    std::vector<Vector3f> centroids_list;
    for (entities::Mesh::VertexFaceIter it_face = mesh.vf_iter(vertex);
         it_face.is_valid(); ++it_face) {
        std::cout << *it_face << " ";

        const auto centroid = face_centroid(mesh, it_face);
        centroids_list.push_back(centroid);
    }
    std::cout << std::endl;

    MatrixXf centroids(centroids_list.size(), 3);
    for (int i = 0; i < centroids_list.size(); ++i) {
        centroids.row(i) = centroids_list[i];
    }

    return centroids;
}

namespace services {

    // IO

    entities::Mesh MeshService::load_mesh(
            const std::string &filename
    ) const {
        spdlog::info("Loading mesh from file {}", filename);
        return this->mesh_dao.load_mesh_from_file(filename);
    }

    bool MeshService::is_trimesh(
            const entities::Mesh &mesh
    ) const {
        spdlog::info("Checking if mesh is a triangle mesh");

        int n_triangles = 0;
        int n_non_triangles = 0;
        for (auto it_f = mesh.faces_begin(); it_f != mesh.faces_end(); ++it_f) {
            if (mesh.valence(*it_f) == 3) {
                n_triangles++;
            } else {
                n_non_triangles++;
            }
        }

        spdlog::debug("Finished checking if the mesh is a triangle mesh with Triangles={}, Non Triangle={}",
                      n_triangles, n_non_triangles);

        return n_non_triangles == 0;
    }

    entities::UnboundModel MeshService::load_unbound_model_from_file(
            const std::string &filename
    ) const {
        spdlog::info("Loading SDFn from file {}", filename);

        return this->mesh_dao.load_unbound_model(filename);
    }

    void MeshService::save_mesh(
            const std::string &filename,
            const Parametrizer &field
    ) const {
        spdlog::info("Saving mesh to file {}", filename);

        std::ofstream os(filename);
        for (int i = 0; i < field.m_positions_compact.size(); ++i) {
            auto t = field.m_positions_compact[i] * field.m_normalize_scale + field.m_normalize_offset;
            os << "v " << t[0] << " " << t[1] << " " << t[2] << "\n";
        }
        for (auto &i: field.m_faces_compact) {
            os << "f " << i[0] + 1 << " " << i[1] + 1 << " "
               << i[2] + 1 << " " << i[3] + 1 << "\n";
        }
        os.close();
    }

    void MeshService::save_mesh(
            const std::string &filename,
            const entities::Mesh &mesh
    ) const {
        spdlog::info("Saving mesh to file {}", filename);

        mesh_dao.save_mesh_to_file(filename, mesh);
    }

    // Quadriflow

    void MeshService::set_boundary_constraints(
            Hierarchy &hierarchy
    ) {
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

    std::map<int, int> MeshService::find_orientation_singularities(
            Hierarchy &hierarchy
    ) {
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
    ) {
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

    // Conversions

    entities::Mesh MeshService::to_trimesh(
            entities::Mesh &mesh
    ) const {
        spdlog::info("Converting mesh to triangle mesh vertices={} faces={}", mesh.n_vertices(), mesh.n_faces());

        spdlog::stopwatch watch;

        mesh.request_face_status();
        for (entities::Mesh::FaceIter it_f = mesh.faces_begin(); it_f != mesh.faces_end(); ++it_f) {
            if (mesh.valence(*it_f) == 3) continue;

            if (mesh.valence(*it_f) == 4) { // quad
                entities::Mesh::CFVIter fv_it = mesh.cfv_iter(*it_f);
                auto v0 = *fv_it;
                auto v1 = *(++fv_it);
                auto v2 = *(++fv_it);
                auto v3 = *(++fv_it);

                mesh.delete_face(*it_f, false);
                mesh.add_face(v0, v1, v2);
                mesh.add_face(v0, v2, v3);
            } else {
                spdlog::warn("Encountered non-quad face which is not yet supported in the triangle mesh conversion");
            }
        }

        mesh.garbage_collection();
        spdlog::debug("Finished converting mesh to triangle mesh vertices={} faces={} ({:.3}s)",
                      mesh.n_vertices(), mesh.n_faces(), watch);

        return mesh;
    }

    // Fields

    Vector3f MeshService::create_laplacian_angle_field(
            const entities::SDFn &sdfn,
            entities::Mesh &mesh
    ) const {

        VectorXf divergences(mesh.n_vertices());

//        tbb::parallel_for(size_t(0), mesh.n_vertices(), [&](size_t idx) {
        for (int idx = 0; idx < mesh.n_vertices(); ++idx) {
            entities::Mesh::VertexHandle it_vertex(idx);

            const auto centroids = face_centroids_ring(mesh, it_vertex);
            std::cout << "centroids:\n" << centroids << std::endl;

            const auto face_normals = sdfn::normal_of(sdfn, centroids);
            std::cout << "face_normals:\n" << face_normals << std::endl;

            MatrixXf angles = face_normals * face_normals.transpose();
            angles = clip(angles, -1.f, 1.f).array().acos() * (180.f / M_PI);
            std::cout << "angles:\n" << angles << std::endl;

            divergences[idx] = frobenius_norm_off_diagonal(angles);
        }

//        });

        return divergences;
    }

    // Smoothing

    entities::Mesh MeshService::smoothing_surface_snapping(
            const entities::SDFn &sdfn,
            entities::Mesh &mesh,
            const int iterations,
            const float rate
    ) const {
        spdlog::info("Smoothing mesh with SDFn gradient. iterations={}, rate={}", iterations, rate);
        spdlog::stopwatch watch;

        auto vertices = MatrixXf(mesh.n_vertices(), 3);
        for (int i = 0; i < mesh.n_vertices(); ++i) {
            auto point = mesh.point(entities::Mesh::VertexHandle(i));
            vertices.row(i) = Vector3f(point[0], point[1], point[2]);
        }

        const float error_before = sdfn(vertices).cwiseAbs().sum();

        for (int i = 0; i < iterations; ++i) {
            spdlog::debug("Smoothing iteration {}/{}", i, iterations);

            const auto gradients = sdfn::gradient_of(sdfn, vertices);
            vertices -= sdfn(vertices) * gradients * rate;
        }

        const float error_after = sdfn(vertices).cwiseAbs().sum();

        for (int i = 0; i < mesh.n_vertices(); ++i) {
            mesh.set_point(
                    entities::Mesh::VertexHandle(i),
                    entities::Mesh::Point(
                            vertices(i, 0),
                            vertices(i, 1),
                            vertices(i, 2)
                    )
            );
        }

        spdlog::debug("Finished smoothing mesh with SDFn gradient distances before={} after={} ({:.3}s)",
                      error_before, error_after, watch);

#ifdef DEV_DEBUG
        OpenMesh::IO::write_mesh(mesh, "../tests/out/smoothing_surface_snapping.ply");
#endif

        return mesh;
    }

    entities::Mesh MeshService::smoothing_edge_snapping(
            const entities::SDFn &sdfn,
            entities::Mesh &mesh,
            const int iterations,
            const float threshold_angle,
            const float max_error
    ) const {
        spdlog::info("Smoothing mesh by snapping vertices to edges.");
        spdlog::stopwatch watch;

        for (int iteration = 0; iteration < iterations; ++iteration) {
            const auto field = create_laplacian_angle_field(sdfn, mesh);

            MatrixXf vertices_smoothed(mesh.n_vertices(), 3);

            int i = 0;
            for (auto it_v = mesh.vertices_begin(); it_v != mesh.vertices_end(); ++it_v) {
                const auto point = mesh.point(*it_v);
                const auto vertex = Vector3f(point[0], point[1], point[2]);
                vertices_smoothed.row(i) = vertex;


                if (field[i] > threshold_angle) {
                    // Optimize: Neighborhood is computed twice
                    const auto centroids = face_centroids_ring(mesh, *it_v);
                    const auto face_normals = sdfn::normal_of(sdfn, centroids);
                    const auto vertex_new = transformations::intersect_planes(centroids, face_normals);

                    if (vertex_new[3] < max_error) {
                        // TODO: Check if new location is outside neighborhood
                        vertices_smoothed.row(i) = vertex_new.head(3);
                    }
                }

                i++;
            }

            // FIXME: Remove once eigen is used as storage
            for (int i = 0; i < mesh.n_vertices(); ++i) {
                mesh.set_point(
                        entities::Mesh::VertexHandle(i),
                        entities::Mesh::Point(
                                vertices_smoothed(i, 0),
                                vertices_smoothed(i, 1),
                                vertices_smoothed(i, 2)
                        )
                );
            }
        }

        spdlog::debug("Finished smoothing mesh by snapping vertices to edges ({:.3}s)", watch);

#ifdef DEV_DEBUG
        OpenMesh::IO::write_mesh(mesh, "../tests/out/smoothing_edge_snapping.ply");
#endif

        return mesh;
    }

    // Top level

    entities::Mesh MeshService::mesh(
            const entities::SDFn &sdfn,
            const int resolution,
            const AlignedBox3f &bounds
    ) const {
        spdlog::info("Meshing SDFn via surface nets with resolution {}", resolution);

        surfacenets::SurfaceNetsMeshStrategy strategy;
        return strategy.mesh(sdfn, resolution, bounds);
    }

    Parametrizer MeshService::remesh(
            const entities::Mesh &mesh,
            const int face_count,
            const bool preserve_edges,
            const bool preserve_boundaries,
            const bool use_adaptive_meshing
    ) const {
        assert(is_trimesh(mesh));

        Parametrizer field;
        spdlog::info("Re-meshing mesh with {} target faces", face_count);

        spdlog::stopwatch watch;

        spdlog::info("Initializing parameters");
        spdlog::debug("Re-meshing mesh with vertices={}, faces={}", mesh.n_vertices(), mesh.n_faces());

        adapters::initialize_parameterizer(field, mesh);
        field.initialize_parameterizer(
                preserve_boundaries,
                preserve_edges,
                face_count,
                use_adaptive_meshing
        );

        spdlog::debug("Finished initializing parameters ({:.3}s)", watch);

        if (preserve_boundaries) {
            set_boundary_constraints(field.m_hierarchy);
        }

        watch.reset();
        spdlog::info("Solving orientation field");

        Optimizer::optimize_orientations(field.m_hierarchy);
        find_orientation_singularities(field.m_hierarchy);

        spdlog::debug("Finished solving orientation field ({:.3}s)", watch);


        if (use_adaptive_meshing) {
            watch.reset();
            spdlog::info("Analyzing mesh for adaptive scaling");

            const auto [faces_slope, faces_orientation] = estimate_slope(
                    field.m_hierarchy,
                    field.m_triangle_space,
                    field.m_faces_normals
            );
            field.m_faces_slope = faces_slope;
            field.m_faces_orientation = faces_orientation;

            spdlog::info("Finished analyzing mesh for adaptive scaling ({:.3}s)", watch);
        }


        watch.reset();
        spdlog::info("Solving field for adaptive scale");

        Optimizer::optimize_scale(field.m_hierarchy, field.m_rho, use_adaptive_meshing);

        spdlog::debug("Finished solving field for adaptive scale ({:.3}s)", watch);


        watch.reset();
        spdlog::info("Solving for position field");

        Optimizer::optimize_positions(field.m_hierarchy);
        const auto [singularity_position, singularity_rank, singularity_index] = find_position_singularities(
                field.m_hierarchy,
                true
        );
        field.m_singularity_position = singularity_position;
        field.m_singularity_rank = singularity_rank;
        field.m_singularity_index = singularity_index;

        spdlog::debug("Finished solving for position field ({:.3}s)", watch);

        watch.reset();
        spdlog::info("Solving index map");

        field.compute_index_map(field.m_hierarchy);

        spdlog::debug("Finished solving index map ({:.3}s)", watch);

        return field;
    }
}