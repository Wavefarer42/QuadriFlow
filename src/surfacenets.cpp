#include <Eigen/Dense>
#include <OpenMesh/Tools/Utils/MeshCheckerT.hh>

#include "surfacenets.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <OpenMesh/Core/IO/MeshIO.hh>

using namespace Eigen;

namespace surfacenets {

    AlignedBox3f SurfaceNetsMeshStrategy::estimate_bounding_box(
            entities::SDFn sdfn,
            int resolution
    ) const {
        spdlog::debug("Estimating bounding box");
        spdlog::stopwatch watch;

        float step_size = 0.5f;

        VectorXf lower_bound = VectorXf::Constant(3, std::numeric_limits<float>::max());
        VectorXf upper_bound = VectorXf::Constant(3, std::numeric_limits<float>::lowest());

        const int total_points = (2 * resolution + 1) * (2 * resolution + 1) * (2 * resolution + 1);

        MatrixXf domain(total_points, 3);
        int index = 0;
        for (int x = -resolution; x <= resolution; ++x) {
            for (int y = -resolution; y <= resolution; ++y) {
                for (int z = -resolution; z <= resolution; ++z) {
                    domain.row(index++) << x * step_size, y * step_size, z * step_size;
                }
            }
        }

        const VectorXf distances = sdfn(domain);

        for (int i = 0; i < total_points; ++i) {
            if (distances[i] <= 0.0f) {
                lower_bound = lower_bound.array().min(domain.row(i).transpose().array());
                upper_bound = upper_bound.array().max(domain.row(i).transpose().array());
            }
        }

        AlignedBox3f bounds(lower_bound, upper_bound);
        spdlog::debug("Estimated bounding box: min = ({}, {}, {}), max = ({}, {}, {}), elapsed = {}",
                      bounds.min().x(), bounds.min().y(), bounds.min().z(),
                      bounds.max().x(), bounds.max().y(), bounds.max().z(), watch);
        return bounds;
    }

    VectorXi SurfaceNetsMeshStrategy::create_face(
            const VectorXi &face_indices,
            const MatrixXf &face_vertices,
            bool is_negative_face
    ) const {
        VectorXi face(4);

        if (is_negative_face) {
            face << face_indices[1],
                    face_indices[3],
                    face_indices[2],
                    face_indices[0];
        } else {
            face << face_indices[0],
                    face_indices[2],
                    face_indices[3],
                    face_indices[1];
        }



//        VectorXi face(4);
//        if (distance14 < distance23) {
//            if (is_negative_face) {
//                face << face_indices[0],
//                        face_indices[3],
//                        face_indices[1],
//                        face_indices[2];
//            } else {
//                face << face_indices[0],
//                        face_indices[1],
//                        face_indices[3],
//                        face_indices[2];
//            }
//        } else if (is_negative_face) {
//            face << face_indices[1],
//                    face_indices[2],
//                    face_indices[3],
//                    face_indices[0];
//        } else {
//            face << face_indices[1],
//                    face_indices[3],
//                    face_indices[2],
//                    face_indices[0];
//        }

        return face;
    }

    Vector3f SurfaceNetsMeshStrategy::estimate_centroid(const VectorXf &distances_corners) const {

        int n = 0;
        Vector3f sum = Vector3f::Zero();
        for (int i = 0; i < CUBE_EDGES.rows(); ++i) {
            const auto idx_edge_a = CUBE_EDGES(i, 0);
            const auto idx_edge_b = CUBE_EDGES(i, 1);
            const auto distance_a = distances_corners(idx_edge_a);
            const auto distance_b = distances_corners(idx_edge_b);

            if ((distance_a < 0.0f) != (distance_b < 0.0f)) {
                const auto interpolation1 = distance_a / (distance_a - distance_b + 1e-10f);
                const auto interpolation2 = 1.0f - interpolation1;

                const Vector3f interpolation = interpolation2 * CUBE_CORNERS.row(idx_edge_a).cast<float>()
                                               + interpolation1 * CUBE_CORNERS.row(idx_edge_b).cast<float>();
                sum += interpolation;
                n++;
            }
        }

        if (n > 0) {
            return sum / static_cast<float>(n);
        } else {
            return Vector3f::Zero();
        }
    }

    entities::Mesh SurfaceNetsMeshStrategy::mesh(
            const entities::SDFn &sdfn,
            const entities::Shape shape
    ) const {
        spdlog::stopwatch watch, watch_total;

        AlignedBox3f bounds_ = shape.bounds;
        if (bounds_.volume() == 0) {
            bounds_ = estimate_bounding_box(sdfn, shape.resolution);
        }

        // Create domains
        watch.reset();
        spdlog::debug("Creating sample domain");

        // [0:3]: Grid coordinate (x, y, z)
        // [3]: Mapping from grid coordinate to linear surface array index see vertices
        MatrixXi indices(shape.size(), 4);
        int index = 0;
        for (int z = 0; z < shape.resolution + 1; ++z) {
            for (int y = 0; y < shape.resolution + 1; ++y) {
                for (int x = 0; x < shape.resolution + 1; ++x) {
                    indices.row(index++) << x, y, z, -1;
                }
            }
        }

        const auto linearize = [shape](Vector3i idx_nd) {
            return idx_nd.x()
                   + idx_nd.y() * (shape.resolution + 1)
                   + idx_nd.z() * (shape.resolution + 1) * (shape.resolution + 1);
        };

        MatrixXf domain = indices.block(0, 0, indices.rows(), 3).cast<float>();
        domain /= static_cast<float>(shape.resolution);
        domain.array().rowwise() *= bounds_.sizes().transpose().array();
        domain.array().rowwise() += bounds_.min().transpose().array();

        spdlog::debug("Finished creating sample domain ({:.3}s)", watch);

        watch.reset();
        spdlog::debug("Sampling signed distance field");

        const MatrixXf sdf = sdfn(domain);

        spdlog::debug("Finished sampling signed distance field ({:.3}s)", watch);


        // Create vertices
        watch.reset();
        spdlog::debug("Finding vertex positions");

        std::vector<Vector3f> vertices;
        for (int i = 0; i < indices.rows(); ++i) {
            if (indices(i, 0) >= shape.resolution
                || indices(i, 1) >= shape.resolution
                || indices(i, 2) >= shape.resolution) {
                continue;
            }

            const Vector3i corner = indices.row(i).head<3>();

            // Get the field at the sample grid corner points
            VectorXf distances_corners(CUBE_CORNERS.rows());
            for (int j = 0; j < CUBE_CORNERS.rows(); ++j) {
                const Vector3i idx_nd = corner + CUBE_CORNERS.row(j).transpose();
                const int idx_flat = linearize(idx_nd);
                distances_corners(j) = sdf(idx_flat);
            }

            auto num_negatives = (distances_corners.array() < 0.0f).count();
            if (num_negatives != 0 && num_negatives != 8) {
                const auto centroid = estimate_centroid(distances_corners);
                vertices.emplace_back(corner.cast<float>() + centroid);
                indices(i, 3) = static_cast<int>(vertices.size()) - 1;
            }
        }

        const auto on_surface = [](float d1, float d2) {
            return (d1 < 0 && 0 < d2) || (d2 < 0 && 0 < d1);
        };
        const auto is_negative_face = [](float d1, float d2) {
            return d1 < 0 && 0 < d2;
        };

        const auto fn_face_vertices = [linearize, vertices, indices](
                const Vector3i &vidx,
                const Vector3i &axis_b,
                const Vector3i &axis_c
        ) {
            Vector4i face_indices(4);
            face_indices[0] = indices(linearize(vidx), 3);
            face_indices[1] = indices(linearize(vidx - axis_b), 3);
            face_indices[2] = indices(linearize(vidx - axis_c), 3);
            face_indices[3] = indices(linearize(vidx - axis_b - axis_c), 3);

            MatrixXf face_vertices(4, 3);
            for (int i = 0; i < face_vertices.rows(); ++i) {
                assert(face_indices[i] >= 0); // Index mapping points to invalid vertex
                face_vertices.row(i) = vertices[face_indices[i]];
            }

            return std::make_pair(face_indices, face_vertices);
        };

        spdlog::debug("Finished finding vertex positions ({:.3}s)", watch);

        // Create faces
        watch.reset();
        spdlog::debug("Creating faces");

        std::vector<VectorXi> faces;
        for (int i = 0; i < indices.rows(); ++i) {
            if (indices(i, 3) == -1) continue;

            Vector3i idx_vertex_1 = indices.row(i).head<3>();

            if (idx_vertex_1.x() != shape.resolution
                && idx_vertex_1.y() != 0
                && idx_vertex_1.z() != 0) {
                const Vector3i idx_vertex_2 = idx_vertex_1 + AXIS_X;
                const auto d1 = sdf(linearize(idx_vertex_1));
                const auto d2 = sdf(linearize(idx_vertex_2));

                if (on_surface(d1, d2)) {
                    const auto [face_indices, face_vertices] = fn_face_vertices(idx_vertex_1, AXIS_Y, AXIS_Z);
                    const auto face = create_face(face_indices, face_vertices, is_negative_face(d1, d2));
                    faces.emplace_back(face);
                }
            }

            if (idx_vertex_1.x() != 0
                && idx_vertex_1.y() != shape.resolution
                && idx_vertex_1.z() != 0) {
                const Vector3i idx_vertex_2 = idx_vertex_1 + AXIS_Y;
                const auto d1 = sdf(linearize(idx_vertex_1));
                const auto d2 = sdf(linearize(idx_vertex_2));

                if (on_surface(d1, d2)) {
                    const auto [face_indices, face_vertices] = fn_face_vertices(idx_vertex_1, AXIS_Z, AXIS_X);
                    const auto face = create_face(face_indices, face_vertices, is_negative_face(d1, d2));
                    faces.emplace_back(face);
                }
            }

            if (idx_vertex_1.x() != 0
                && idx_vertex_1.y() != 0
                && idx_vertex_1.z() != shape.resolution) {
                const Vector3i idx_vertex_2 = idx_vertex_1 + AXIS_Z;
                const auto d1 = sdf(linearize(idx_vertex_1));
                const auto d2 = sdf(linearize(idx_vertex_2));

                if (on_surface(d1, d2)) {
                    const auto [face_indices, face_vertices] = fn_face_vertices(idx_vertex_1, AXIS_X, AXIS_Y);
                    const auto face = create_face(face_indices, face_vertices, is_negative_face(d1, d2));
                    faces.emplace_back(face);
                }
            }
        }

        spdlog::debug("Finished creating faces ({:.3}s)", watch);


        // Create final mesh;
        watch.reset();
        spdlog::debug("Creating mesh data structure");

        entities::Mesh mesh;
        std::vector<entities::Mesh::VertexHandle> vertex_handles;
        for (const auto &vertex: vertices) {
            const auto handle = mesh.add_vertex(entities::Mesh::Point(vertex.x(), vertex.y(), vertex.z()));
            vertex_handles.emplace_back(handle);
        }

        for (const auto &face: faces) {
            mesh.add_face(vertex_handles[face[0]],
                          vertex_handles[face[1]],
                          vertex_handles[face[2]],
                          vertex_handles[face[3]]);
        }

        assert(OpenMesh::IO::write_mesh(mesh, "/Users/hannes/Projects/QuadriFlow/tests/resources/sphere.ply"));


        spdlog::debug("Finished creating mesh data structure ({:.3}s)", watch);

        OpenMesh::Utils::MeshCheckerT<entities::Mesh> checker(mesh);
        if (checker.check()) {
            spdlog::debug("Mesh is valid");
        } else {
            spdlog::error("Mesh is invalid");
        }

        return mesh;
    }

}