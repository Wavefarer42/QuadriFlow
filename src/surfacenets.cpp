#include <Eigen/Dense>
#include <OpenMesh/Tools/Utils/MeshCheckerT.hh>

#include "surfacenets.h"
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"
#include <OpenMesh/Core/IO/MeshIO.hh>

using namespace Eigen;

namespace surfacenets {

    bool is_on_surface(
            const float distance1,
            const float distance2
    ) {
        return (distance1 < 0 && 0 < distance2) || (distance2 < 0 && 0 < distance1);
    }

    bool is_negative_face(
            const float distance1,
            const float distance2
    ) {
        return distance1 < 0 && 0 < distance2;
    }

    AlignedBox3f SurfaceNetsMeshStrategy::estimate_bounding_box(
            const entities::SDFn &sdfn,
            int resolution
    ) {
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
                    domain.row(index++) << static_cast<float>(x) * step_size,
                            static_cast<float>(y) * step_size,
                            static_cast<float>(z) * step_size;
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

    VectorXi create_face(
            const VectorXi &face_indices,
            bool is_negative_face
    ) {
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

        return face;
    }

    Vector3f SurfaceNetsMeshStrategy::estimate_centroid(
            const VectorXf &distances_corners
    ) const {

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

    MatrixXi create_indices(int resolution) {

        spdlog::debug("Creating domain indices");
        spdlog::stopwatch watch;

        // [0:3]: Grid coordinate (x, y, z)
        // [3]: Mapping from grid coordinate to linear surface array index see vertices
        const auto grid_max = resolution + 1;
        MatrixXi indices(grid_max * grid_max * grid_max, 4);
        int index = 0;
        for (int z = 0; z < grid_max; ++z) {
            for (int y = 0; y < grid_max; ++y) {
                for (int x = 0; x < grid_max; ++x) {
                    indices.row(index++) << x, y, z, -1;
                }
            }
        }

        spdlog::debug("- Created domain indices ({:.3}s)", watch);

        return indices;
    }

    MatrixXf scale_to_domain(
            const MatrixXi &indices,
            const AlignedBox3f &bounds,
            const int resolution
    ) {

        spdlog::debug("Scaling domain to bounding box and resolution");
        spdlog::stopwatch watch;

        MatrixXf domain = indices.block(0, 0, indices.rows(), 3).cast<float>();
        domain /= static_cast<float>(resolution);
        domain.array().rowwise() *= bounds.sizes().transpose().array();
        domain.array().rowwise() += bounds.min().transpose().array();

        spdlog::debug("- Scaled domain to bounding box and resolution ({:.3}s)", watch);

        return domain;
    }

    Vector4i gather_face_indices(
            const Vector3i &vidx,
            const Vector3i &axis_b,
            const Vector3i &axis_c,
            const MatrixXi &indices,
            const NdToFlatIndexer &linearize
    ) {
        Vector4i face_indices(4);
        face_indices[0] = indices(linearize(vidx), 3);
        face_indices[1] = indices(linearize(vidx - axis_b), 3);
        face_indices[2] = indices(linearize(vidx - axis_c), 3);
        face_indices[3] = indices(linearize(vidx - axis_b - axis_c), 3);

        assert(face_indices.minCoeff() >= 0);

        return face_indices;
    }

    MatrixXf sample_sdf(
            const entities::SDFn &sdfn,
            const MatrixXf &domain
    ) {

        spdlog::debug("Sampling signed distance field");
        spdlog::stopwatch watch;

        const MatrixXf sdf = sdfn(domain);

        spdlog::debug("- Sampled signed distance field ({:.3}s)", watch);

        return sdf;
    }

    std::vector<Vector3f> SurfaceNetsMeshStrategy::create_vertices(
            MatrixXi &indices,
            const MatrixXf &sdf,
            const int resolution,
            const NdToFlatIndexer &linearize
    ) const {
        spdlog::debug("Creating surface vertices");
        spdlog::stopwatch watch;

        std::vector<Vector3f> vertices;
        for (int i = 0; i < indices.rows(); ++i) {
            if (indices(i, 0) >= resolution
                || indices(i, 1) >= resolution
                || indices(i, 2) >= resolution) {
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

        spdlog::debug("- Created surface vertices ({:.3}s)", watch);

        return vertices;
    }

    std::vector<VectorXi> SurfaceNetsMeshStrategy::create_faces(
            const MatrixXi &indices,
            const MatrixXf &sdf,
            const int resolution,
            const NdToFlatIndexer &linearize
    ) const {

        spdlog::debug("Creating faces");
        spdlog::stopwatch watch;

        std::vector<VectorXi> faces;
        for (int i = 0; i < indices.rows(); ++i) {
            if (indices(i, 3) == -1) continue;

            Vector3i idx_vertex_1 = indices.row(i).head<3>();

            if (idx_vertex_1.x() != resolution
                && idx_vertex_1.y() != 0
                && idx_vertex_1.z() != 0) {
                const Vector3i idx_vertex_2 = idx_vertex_1 + AXIS_X;
                const auto d1 = sdf(linearize(idx_vertex_1));
                const auto d2 = sdf(linearize(idx_vertex_2));

                if (is_on_surface(d1, d2)) {
                    const auto face_indices = gather_face_indices(idx_vertex_1, AXIS_Y, AXIS_Z, indices, linearize);
                    const auto face = create_face(face_indices, is_negative_face(d1, d2));
                    faces.emplace_back(face);
                }
            }

            if (idx_vertex_1.x() != 0
                && idx_vertex_1.y() != resolution
                && idx_vertex_1.z() != 0) {
                const Vector3i idx_vertex_2 = idx_vertex_1 + AXIS_Y;
                const auto d1 = sdf(linearize(idx_vertex_1));
                const auto d2 = sdf(linearize(idx_vertex_2));

                if (is_on_surface(d1, d2)) {
                    const auto face_indices = gather_face_indices(idx_vertex_1, AXIS_Z, AXIS_X, indices, linearize);
                    const auto face = create_face(face_indices, is_negative_face(d1, d2));
                    faces.emplace_back(face);
                }
            }

            if (idx_vertex_1.x() != 0
                && idx_vertex_1.y() != 0
                && idx_vertex_1.z() != resolution) {
                const Vector3i idx_vertex_2 = idx_vertex_1 + AXIS_Z;
                const auto d1 = sdf(linearize(idx_vertex_1));
                const auto d2 = sdf(linearize(idx_vertex_2));

                if (is_on_surface(d1, d2)) {
                    const auto face_indices = gather_face_indices(idx_vertex_1, AXIS_X, AXIS_Y, indices, linearize);
                    const auto face = create_face(face_indices, is_negative_face(d1, d2));
                    faces.emplace_back(face);
                }
            }
        }

        spdlog::debug("- Finished creating faces ({:.3}s)", watch);

        return faces;
    }

    entities::Mesh finalize_mesh(
            const std::vector<Vector3f> &vertices,
            const std::vector<VectorXi> &faces
    ) {

        spdlog::debug("Creating mesh data structure");
        spdlog::stopwatch watch;

        entities::Mesh mesh;
        std::vector<entities::Mesh::VertexHandle> handles;
        for (const auto &vertex: vertices) {
            const auto handle = mesh.add_vertex(entities::Mesh::Point(vertex.x(), vertex.y(), vertex.z()));
            handles.emplace_back(handle);
        }

        for (const auto &face: faces) {
            mesh.add_face(handles[face[0]], handles[face[1]], handles[face[2]], handles[face[3]]);
        }

        OpenMesh::Utils::MeshCheckerT<entities::Mesh> checker(mesh);
        if (checker.check()) {
            spdlog::debug("Mesh is valid");
        } else {
            spdlog::error("Mesh is invalid");
        }

        spdlog::debug("- Created mesh data structure ({:.3}s)", watch);

        return mesh;
    }

    entities::Mesh SurfaceNetsMeshStrategy::mesh(
            const entities::SDFn &sdfn,
            const AlignedBox3f &bounds,
            const int resolution
    ) const {
        spdlog::stopwatch watch, watch_total;

        AlignedBox3f bounds_ = bounds;
        if (bounds_.volume() == 0) {
            bounds_ = estimate_bounding_box(sdfn, resolution);
        }

        MatrixXi indices = create_indices(resolution);

        const auto linearize = [resolution](Vector3i idx_nd) {
            return idx_nd.x()
                   + idx_nd.y() * (resolution + 1)
                   + idx_nd.z() * (resolution + 1) * (resolution + 1);
        };

        const MatrixXf domain = scale_to_domain(indices, bounds_, resolution);
        const MatrixXf sdf = sample_sdf(sdfn, domain);
        const std::vector<Vector3f> vertices = create_vertices(indices, sdf, resolution, linearize);
        const std::vector<VectorXi> faces = create_faces(indices, sdf, resolution, linearize);
        const entities::Mesh mesh = finalize_mesh(vertices, faces);

        return mesh;
    }

}