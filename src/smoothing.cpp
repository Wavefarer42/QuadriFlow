#include "smoothing.h"

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/IO/Options.hh>

#include "mathext.h"
#include "sdfn.h"
#include "cmake-build-release/_deps/spdlog-src/include/spdlog/stopwatch.h"
#include "spdlog/spdlog.h"

namespace smoothing {
    // Fields

    VectorXf laplacian_angular_field(
        const entities::SDFn &sdfn,
        entities::Mesh &mesh
    ) {
        VectorXf divergences(mesh.n_vertices());

        tbb::parallel_for(size_t(0), mesh.n_vertices(), [&](size_t idx) {
            entities::Mesh::VertexHandle it_vertex(idx);

            const auto centroids = mathext::face_centroids_ring(mesh, it_vertex);
            const auto face_normals = sdfn::normal_of(sdfn, centroids);

            MatrixXf angles = face_normals * face_normals.transpose();
            angles = mathext::clip(angles, -1.f, 1.f).array().acos() * (180.f / M_PI);

            divergences[idx] = angles.maxCoeff();
        });

#ifdef DEV_DEBUG
        entities::Mesh mesh_face_normals;

        mesh_face_normals.request_vertex_normals();
        for (auto it = mesh.faces_begin(); it != mesh.faces_end(); ++it) {
            const auto centroid = mathext::face_centroid(mesh, *it);

            const OpenMesh::VertexHandle vh = mesh_face_normals.add_vertex(
                entities::Mesh::Point(centroid[0], centroid[1], centroid[2])
            );
            Vector3f normal = sdfn::normal_of(sdfn, centroid.matrix().transpose()).row(0);
            mesh_face_normals.set_normal(vh, entities::Mesh::Normal(normal[0], normal[1], normal[2]));
        }

        assert(mesh_face_normals.has_vertex_normals());

        OpenMesh::IO::Options options;
        options += OpenMesh::IO::Options::VertexNormal;
        OpenMesh::IO::write_mesh(mesh_face_normals, "../tests/out/stage/7-laplacian-angle-field.ply", options);

#endif

        return divergences;
    }

    // Smoothing

    entities::Mesh sdfn_projection(
        const entities::SDFn &sdfn,
        entities::Mesh &mesh,
        const int iterations,
        const float rate
    ) {
        spdlog::info("Smoothing mesh with SDFn gradient. iterations={}, rate={}", iterations, rate);
        spdlog::stopwatch watch;

        auto vertices = MatrixXf(mesh.n_vertices(), 3);
        for (int i = 0; i < mesh.n_vertices(); ++i) {
            auto point = mesh.point(entities::Mesh::VertexHandle(i));
            vertices.row(i) = Vector3f(point[0], point[1], point[2]);
        }

        const float error_before = sdfn(vertices).cwiseAbs().sum();

        for (int i = 0; i < iterations; ++i) {
            spdlog::debug("Surface Smoothing iteration {}/{}", i, iterations);

            const VectorXf scale = sdfn(vertices) * rate;
            const MatrixXf direction = sdfn::normal_of(sdfn, vertices);
            const MatrixXf update = direction.array().colwise() * scale.array();

            vertices -= update;
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
        OpenMesh::IO::write_mesh(mesh, "../tests/out/stage/smoothing_surface_snapping.ply");
#endif

        return mesh;
    }

    entities::Mesh edge_snapping(
        const entities::SDFn &sdfn,
        entities::Mesh &mesh,
        const int iterations,
        const float threshold_angle,
        const float max_error
    ) {
        spdlog::info("Smoothing mesh by snapping vertices to edges.");
        spdlog::stopwatch watch;


        for (int iteration = 0; iteration < iterations; ++iteration) {
            const auto field = laplacian_angular_field(sdfn, mesh);

            MatrixXf vertices_smoothed(mesh.n_vertices(), 3);

            tbb::parallel_for(static_cast<size_t>(0), mesh.n_vertices(), [&](size_t idx) {
                const auto it_v = mesh.vertex_handle(idx);
                const auto point = mesh.point(it_v);
                const auto vertex = Vector3f(point[0], point[1], point[2]);
                vertices_smoothed.row(idx) = vertex;

                if (field[idx] > threshold_angle) {
                    const auto centroids = mathext::face_centroids_ring(mesh, it_v);
                    const auto face_normals = sdfn::normal_of(sdfn, centroids);
                    const auto xerr = mathext::intersect_planes(centroids, face_normals);

                    if (xerr[3] < max_error) {
                        // Move a bit into direction of the new position
                        // FIXME Point might move outside the neighborhood of the face
                        vertices_smoothed.row(idx) = xerr.head<3>();
                    }
                }
            });

            // Update mesh
            tbb::parallel_for(static_cast<size_t>(0), mesh.n_vertices(), [&](size_t idx) {
                mesh.set_point(
                    mesh.vertex_handle(idx),
                    entities::Mesh::Point(
                        vertices_smoothed(idx, 0),
                        vertices_smoothed(idx, 1),
                        vertices_smoothed(idx, 2)
                    )
                );
            });
        }

        spdlog::debug("Finished smoothing mesh by snapping vertices to edges ({:.3}s)", watch);

#ifdef DEV_DEBUG
        OpenMesh::IO::write_mesh(mesh, "../tests/out/stage/smoothing_edge_snapping.ply");
#endif

        return mesh;
    }

    entities::Mesh laplacian_with_sdfn_projection(
        const entities::SDFn &sdfn,
        entities::Mesh &mesh,
        const int iterations,
        const float rate_projection
    ) {
        spdlog::info("Smoothing mesh by laplacian smoothing with SDFn projection.");
        spdlog::stopwatch watch;

        auto vertices = MatrixXf(mesh.n_vertices(), 3);
        for (int i = 0; i < mesh.n_vertices(); ++i) {
            auto point = mesh.point(entities::Mesh::VertexHandle(i));
            vertices.row(i) = Vector3f(point[0], point[1], point[2]);
        }

        for (auto iteration = 0; iteration < iterations; ++iteration) {
            spdlog::trace("Laplacian SDFn projection iteration {}/{}", iteration, iterations);

            for (auto idx_v = 0; idx_v < vertices.rows(); ++idx_v) {
                int support = 0;
                auto laplacian = MatrixXf(1, 3);
                laplacian.setZero();

                for (auto it_vv = mesh.vv_iter(OpenMesh::VertexHandle(idx_v)); it_vv.is_valid(); ++it_vv) {
                    const auto p = mesh.point(*it_vv);
                    laplacian.row(0) += Vector3f(p[0], p[1], p[2]);
                    support++;
                }

                laplacian.row(0).array() /= static_cast<float>(support);

                // Projection
                const VectorXf scale = sdfn(laplacian) * rate_projection;
                const MatrixXf direction = sdfn::normal_of(sdfn, laplacian);
                vertices.row(idx_v) = laplacian - scale * direction;
            }
        }

        // Update mesh
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

        spdlog::debug("Finished smoothing mesh by laplacian SDFn projection ({:.3}s)", watch);
#ifdef DEV_DEBUG
        OpenMesh::IO::write_mesh(mesh, "../tests/out/stage/smoothing-laplacian-SDFn-projection.ply");
#endif

        return mesh;
    }
}
