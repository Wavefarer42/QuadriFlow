#include <OpenMesh/Tools/Utils/MeshCheckerT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <tbb/mutex.h>
#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

#include "surfacenets.h"

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

        spdlog::debug("Finished creating domain indices ({:.3}s)", watch);

#ifdef DEV_DEBUG
        entities::Mesh mesh;
        for (int i = 0; i < indices.rows(); ++i) {
            mesh.add_vertex(entities::Mesh::Point(indices(i, 0), indices(i, 1), indices(i, 2)));
        }
        OpenMesh::IO::write_mesh(mesh, "../tests/out/stage/1-indices.ply");
#endif

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

        spdlog::debug("Finished scaling domain to bounding box and resolution ({:.3}s)", watch);

#ifdef DEV_DEBUG
        entities::Mesh mesh;
        for (int i = 0; i < indices.rows(); ++i) {
            mesh.add_vertex(entities::Mesh::Point(domain(i, 0), domain(i, 1), domain(i, 2)));
        }
        OpenMesh::IO::write_mesh(mesh, "../tests/out/stage/2-domain.ply");
#endif

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

        const VectorXf sdf = sdfn(domain);

        spdlog::debug("Finished sampling signed distance field ({:.3}s)", watch);

#ifdef DEV_DEBUG
        spdlog::debug("SDF contains inside={}, outside={} points.", (sdf.array() < 0.0f).count(),
                      (sdf.array() > 0.0f).count());

        entities::Mesh mesh;
        for (int i = 0; i < domain.rows(); ++i) {
            if (sdf(i) < 0.0f) {
                mesh.add_vertex(entities::Mesh::Point(domain(i, 0), domain(i, 1), domain(i, 2)));
            }
        }
        OpenMesh::IO::write_mesh(mesh, "../tests/out/stage/3-sdf.ply");
#endif

        return sdf;
    }

    std::vector<Vector3f> SurfaceNetsMeshStrategy::create_vertices(
        MatrixXi &indices,
        const MatrixXf &domain,
        const MatrixXf &sdf,
        const int resolution,
        const AlignedBox3f &bounds,
        const NdToFlatIndexer &linearize
    ) const {
        spdlog::debug("Creating surface vertices");
        spdlog::stopwatch watch;

        const Vector3f domain_offset = domain.row(0);
        const Vector3f domain_scale = (domain.row(domain.rows() - 1) - domain.row(0)).array() / resolution;
        std::vector<Vector3f> vertices;

        tbb::mutex mtx;
        tbb::parallel_for(
            tbb::blocked_range<int>(0, static_cast<int>(indices.rows())),
            [&](const tbb::blocked_range<int> &range) {
                for (int i = range.begin(); i < range.end(); ++i) {
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
                        const auto centroid = estimate_centroid(distances_corners) + corner.cast<float>();

                        const Vector3f position =
                                domain_offset + (centroid.array() * domain_scale.array()).matrix();

                        tbb::mutex::scoped_lock lock(mtx);

                        vertices.emplace_back(position);
                        indices(i, 3) = static_cast<int>(vertices.size()) - 1;

                        lock.release();
                    }
                }
            }
        );

        spdlog::debug("Finished creating surface vertices={} ({:.3}s)", vertices.size(), watch);

#ifdef DEV_DEBUG
        entities::Mesh mesh;
        for (int i = 0; i < vertices.size(); ++i) {
            mesh.add_vertex(entities::Mesh::Point(vertices[i].x(), vertices[i].y(), vertices[i].z()));
        }
        OpenMesh::IO::write_mesh(mesh, "../tests/out/stage/4-vertices.ply");
#endif

        return vertices;
    }

    std::vector<Vector4i> SurfaceNetsMeshStrategy::create_faces(
        const MatrixXi &indices,
        const MatrixXf &sdf,
        const int resolution,
        const NdToFlatIndexer &linearize
    ) const {
        spdlog::debug("Creating faces");
        spdlog::stopwatch watch;

        std::vector<Vector4i> faces;
        for (int i = 0; i < indices.rows(); ++i) {
            if (indices(i, 3) == -1) continue;

            Vector3i idx_vertex = indices.row(i).head<3>();

            // Scan along the X-axis
            if (idx_vertex.x() < resolution
                && idx_vertex.y() > 0
                && idx_vertex.z() > 0) {
                const auto d0 = sdf(linearize(idx_vertex));
                const auto d1 = sdf(linearize(idx_vertex + AXIS_X));

                if (is_on_surface(d0, d1)) {
                    const auto face_indices = gather_face_indices(idx_vertex, AXIS_Y, AXIS_Z, indices, linearize);
                    const auto face = create_face(face_indices, is_negative_face(d0, d1));
                    faces.emplace_back(face);
                }
            }

            // Scan along the Y-axis
            if (idx_vertex.x() > 0
                && idx_vertex.y() < resolution
                && idx_vertex.z() > 0) {
                const auto d0 = sdf(linearize(idx_vertex));
                const auto d1 = sdf(linearize(idx_vertex + AXIS_Y));

                if (is_on_surface(d0, d1)) {
                    const auto face_indices = gather_face_indices(idx_vertex, AXIS_Z, AXIS_X, indices, linearize);
                    const auto face = create_face(face_indices, is_negative_face(d0, d1));
                    faces.emplace_back(face);
                }
            }

            // Scan along the Z-axis
            if (idx_vertex.x() > 0
                && idx_vertex.y() > 0
                && idx_vertex.z() < resolution) {
                const auto d0 = sdf(linearize(idx_vertex));
                const auto d1 = sdf(linearize(idx_vertex + AXIS_Z));

                if (is_on_surface(d0, d1)) {
                    const auto face_indices = gather_face_indices(idx_vertex, AXIS_X, AXIS_Y, indices, linearize);
                    const auto face = create_face(face_indices, is_negative_face(d0, d1));
                    faces.emplace_back(face);
                }
            }
        }

        spdlog::debug("Finished creating faces={} ({:.3}s)", faces.size(), watch);

        return faces;
    }

    entities::Mesh finalize_mesh(
        const std::vector<Vector3f> &vertices,
        const std::vector<Vector4i> &faces
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

        spdlog::debug("Finished creating mesh data structure ({:.3}s)", watch);

#ifdef DEV_DEBUG
        OpenMesh::IO::write_mesh(mesh, "../tests/out/stage/5-faces.ply");
#endif

        return mesh;
    }

    entities::Mesh SurfaceNetsMeshStrategy::mesh(
        const entities::SDFn &sdfn,
        const AlignedBox3f &bounds,
        const int resolution
    ) const {
        spdlog::stopwatch watch;

        if (bounds.volume() == 0.0f) {
            throw std::invalid_argument("Bounds must have a volume greater than zero");
        }

        spdlog::info("Creating mesh from SDF with resolution={} and bounds=([{}, {}, {}], [{}, {}, {}])",
                     resolution, bounds.min().x(), bounds.min().y(), bounds.min().z(), bounds.max().x(),
                     bounds.max().y(), bounds.max().z());

        MatrixXi indices = create_indices(resolution);

        const auto linearize = [resolution](Vector3i idx_nd) {
            return idx_nd.x()
                   + idx_nd.y() * (resolution + 1)
                   + idx_nd.z() * (resolution + 1) * (resolution + 1);
        };

        const auto domain = scale_to_domain(indices, bounds, resolution);
        const auto sdf = sample_sdf(sdfn, domain);
        const auto vertices = create_vertices(indices, domain, sdf, resolution, bounds, linearize);
        const auto faces = create_faces(indices, sdf, resolution, linearize);
        const auto mesh = finalize_mesh(vertices, faces);

        spdlog::debug("Finished creating mesh ({})", watch);

        return mesh;
    }
}
