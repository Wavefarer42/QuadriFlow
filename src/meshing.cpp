#include <OpenMesh/Tools/Utils/MeshCheckerT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <tbb/mutex.h>
#include <fstream>
#include <vector>
#include <nlohmann/json.hpp>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>

#include "spdlog/spdlog.h"
#include "spdlog/stopwatch.h"

#include "meshing.h"
#include "mathext.h"
#include "quadriflow.h"

using namespace Eigen;


namespace surfacenets {
    using NdToFlatIndexer = std::function<int(const Vector3i &)>;

    static const MatrixXi CUBE_CORNERS = (
        MatrixXi(8, 3) <<
        0, 0, 0,
        1, 0, 0,
        0, 1, 0,
        1, 1, 0,
        0, 0, 1,
        1, 0, 1,
        0, 1, 1,
        1, 1, 1
    ).finished();

    static const Matrix<int, 12, 2> CUBE_EDGES = (Matrix<int, 12, 2>() << 0b000, 0b001,
                                                  0b000, 0b010,
                                                  0b000, 0b100,
                                                  0b001, 0b011,
                                                  0b001, 0b101,
                                                  0b010, 0b011,
                                                  0b010, 0b110,
                                                  0b011, 0b111,
                                                  0b100, 0b101,
                                                  0b100, 0b110,
                                                  0b101, 0b111,
                                                  0b110, 0b111).finished();

    static const Vector3i AXIS_X = {1, 0, 0};
    static const Vector3i AXIS_Y = {0, 1, 0};
    static const Vector3i AXIS_Z = {0, 0, 1};

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

    Vector3f estimate_centroid(
        const VectorXf &distances_corners
    ) {
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

    MatrixXi create_indices(
        int resolution
    ) {
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

        assert(face_indices.minCoeff() >= 0 && "Inconsistency between face indices from the vertex creation.");

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

    std::vector<Vector3f> create_vertices(
        MatrixXi &indices,
        const MatrixXf &domain,
        const MatrixXf &sdf,
        const int resolution,
        const NdToFlatIndexer &linearize
    ) {
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

    std::vector<Vector4i> create_faces(
        const MatrixXi &indices,
        const MatrixXf &sdf,
        const int resolution,
        const NdToFlatIndexer &linearize
    ) {
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

    std::function<int(Vector3i)> indexer_nd_to_linear(
        int resolution
    ) {
        return [resolution](Vector3i idx_nd) {
            return idx_nd.x()
                   + idx_nd.y() * (resolution + 1)
                   + idx_nd.z() * (resolution + 1) * (resolution + 1);
        };
    }

    entities::Mesh mesh(
        const entities::SDFn &sdfn,
        const AlignedBox3f &bounds,
        const int resolution
    ) {
        spdlog::stopwatch watch;

        if (bounds.volume() == 0.0f) {
            throw std::invalid_argument("Bounds must have a volume greater than zero");
        }

        spdlog::info("Creating mesh from SDF with resolution={} and bounds=([{}, {}, {}], [{}, {}, {}])",
                     resolution, bounds.min().x(), bounds.min().y(), bounds.min().z(), bounds.max().x(),
                     bounds.max().y(), bounds.max().z());

        MatrixXi indices = create_indices(resolution);

        const auto linearize = indexer_nd_to_linear(resolution);
        const auto domain = scale_to_domain(indices, bounds, resolution);
        const auto sdf = sample_sdf(sdfn, domain);
        const auto vertices = create_vertices(indices, domain, sdf, resolution, linearize);
        const auto faces = create_faces(indices, sdf, resolution, linearize);
        const auto mesh = finalize_mesh(vertices, faces);

        spdlog::debug("Finished creating mesh ({})", watch);

        return mesh;
    }
}

namespace delaunay {
    typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
    typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
    typedef Tr::Geom_traits GT;
    typedef GT::Sphere_3 Sphere_3;
    typedef GT::Point_3 Point_3;
    typedef GT::FT FT;

    typedef FT (*Function)(Point_3);

    typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;
    typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

    entities::Mesh mesh(
        const entities::SDFn &sdfn,
        const AlignedBox3f &bounds,
        const int resolution
    ) {
        const auto _sdfn = [sdfn](const Point_3 &p) {
            const Vector3f domain = {
                static_cast<float>(p.x()),
                static_cast<float>(p.y()),
                static_cast<float>(p.z())
            };
            const auto distance = sdfn(domain.transpose());
            return distance(0);
        };

        Tr tr;
        C2t3 c2t3(tr);
        Surface_mesh _mesh;

        const auto diameter = (bounds.max() - bounds.min()).norm();
        const auto radius = diameter / 2;
        constexpr auto angle_bound = 20.0;
        const auto radius_bound = diameter / resolution;
        const auto distance_bound = diameter / 1000;

        spdlog::debug(
            "Meshing with Delaunay triangulation diameter={}, radius={}, angle={}, radius_bound={}, distance_bound={}",
            diameter, radius, angle_bound, radius_bound, distance_bound
        );

        const auto surface = Surface_3(_sdfn, Sphere_3(CGAL::ORIGIN, radius * radius));
        const auto criteria = CGAL::Surface_mesh_default_criteria_3<Tr>(angle_bound, radius_bound, distance_bound);

        CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
        CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, _mesh);

        entities::Mesh mesh;
        std::map<Surface_mesh::Vertex_index, entities::Mesh::VertexHandle> mapping;
        for (Surface_mesh::Vertex_index vi: _mesh.vertices()) {
            const auto p = _mesh.point(vi);
            mapping[vi] = mesh.add_vertex(entities::Mesh::Point(p.x(), p.y(), p.z()));
        }

        for (Surface_mesh::Face_index fi: _mesh.faces()) {
            std::vector<entities::Mesh::VertexHandle> face;
            for (Surface_mesh::Vertex_index vi: _mesh.vertices_around_face(_mesh.halfedge(fi))) {
                face.emplace_back(mapping[vi]);
            }
            mesh.add_face(face);
        }

        return mesh;
    }
}

namespace meshing {
    void initialize_parameterizer(quadriflow::Parametrizer &field, entities::Mesh mesh) {
        spdlog::debug("Initializing parameters");

        field.m_vertices = MatrixXd(3, mesh.n_vertices());
        for (auto it_v = mesh.vertices_begin(); it_v != mesh.vertices_end(); ++it_v) {
            const auto point = mesh.point(*it_v);
            field.m_vertices.col(it_v->idx()) = Vector3d(point[0], point[1], point[2]);
        }

        field.m_faces = MatrixXi(3, mesh.n_faces());
        for (auto it_f = mesh.faces_begin(); it_f != mesh.faces_end(); ++it_f) {
            auto fv_it = mesh.cfv_iter(*it_f);
            for (int i = 0; i < 3; ++i) {
                field.m_faces(i, it_f->idx()) = fv_it->idx();
                ++fv_it;
            }
        }
    }

    entities::Mesh from_parametrizer_to_quad_mesh(const quadriflow::Parametrizer &field) {
        spdlog::info("Converting parametrizer to mesh");

        entities::Mesh mesh;

        for (int i = 0; i < field.m_positions_compact.size(); ++i) {
            const Vector3f t = (field.m_positions_compact[i] * field.m_normalize_scale + field.m_normalize_offset).cast<
                float>();
            mesh.add_vertex(entities::Mesh::Point(t[0], t[1], t[2]));
        }

        for (const auto &i: field.m_faces_compact) {
            mesh.add_face({
                entities::Mesh::VertexHandle(i[0]),
                entities::Mesh::VertexHandle(i[1]),
                entities::Mesh::VertexHandle(i[2]),
                entities::Mesh::VertexHandle(i[3])
            });
        }

        return mesh;
    }

    std::map<int, int> find_orientation_singularities(
        quadriflow::Hierarchy &hierarchy
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
                auto value = mathext::compat_orientation_extrinsic_index_4(
                    Q.col(i),
                    normals.col(i),
                    Q.col(j),
                    normals.col(j)
                );
                index += value.second - value.first;
                abs_index += std::abs(value.second - value.first);
            }
            int index_mod = mathext::modulo(index, 4);
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

    std::tuple<std::map<int, Vector2i>, MatrixXi, MatrixXi> find_position_singularities(
        quadriflow::Hierarchy &m_hierarchy,
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
                Vector3d v0 = mathext::rotate90_by(q[0], n[0], i);
                for (int j = 0; j < 4; ++j) {
                    Vector3d v1 = mathext::rotate90_by(q[1], n[1], j);
                    for (int k = 0; k < 4; ++k) {
                        Vector3d v2 = mathext::rotate90_by(q[2], n[2], k);
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
            for (int k = 0; k < 3; ++k) q[k] = mathext::rotate90_by(q[k], n[k], best[k]);

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
                std::pair<Vector2i, Vector2i> value = mathext::compat_position_extrinsic_index_4(
                    v[k], n[k], q[k], o[k], v[kn], n[kn], q[kn], o[kn], scale_x, scale_y, inv_scale_x,
                    inv_scale_y, scale_x_1, scale_y_1, inv_scale_x_1, inv_scale_y_1, nullptr);
                auto diff = value.first - value.second;
                index += diff;
                singularity_index(k * 2, f) = diff[0];
                singularity_index(k * 2 + 1, f) = diff[1];
            }

            if (index != Vector2i::Zero()) {
                singularity_position[f] = mathext::rshift90(index, best[0]);
            }
        }

        return std::make_tuple(singularity_position, singularity_rank, singularity_index);
    }

    bool is_trimesh(
        const entities::Mesh &mesh
    ) {
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

    entities::Mesh mesh_to_trimesh(
        const entities::SDFn &sdfn,
        const AlignedBox3f &bounds,
        const int resolution,
        const std::string &algorithm
    ) {
        spdlog::info("Meshing SDFn to trimesh. algorithm={}, resolution={}", algorithm, resolution);

        if (algorithm == "delaunay") {
            return delaunay::mesh(sdfn, bounds, resolution);
        }

        throw std::invalid_argument("Invalid algorithm");
    }

    /**
     * Returns an irregular quad mesh from an implicit surface using resolution as base sampling rate.
     * @param sdfn signed distance function
     * @param bounds bounding box
     * @param resolution samplesing resolution
     * @param algorithm currently: surface-nets
     * @return The quad mesh
     */
    entities::Mesh mesh_to_quadmesh(
        const entities::SDFn &sdfn,
        const AlignedBox3f &bounds,
        const int resolution,
        const std::string &algorithm
    ) {
        spdlog::info("Meshing SDFn to quadmesh. algorithm={} resolution {}", algorithm, resolution);

        if (algorithm == "surface-nets") {
            return surfacenets::mesh(sdfn, bounds, resolution);
        }

        throw std::invalid_argument("Invalid algorithm");
    }

    /**
     * Re-meshes a given mesh to a triangle mesh (currently only supports quad mesh inputs).
     * @param mesh A quad mesh
     * @return A triangle mesh
     */
    entities::Mesh remesh_to_trimesh(
        entities::Mesh &mesh
    ) {
        spdlog::info("Re-meshing to triangle mesh vertices={} faces={}", mesh.n_vertices(), mesh.n_faces());

        spdlog::stopwatch watch;

        mesh.request_face_status();
        for (entities::Mesh::FaceIter it_f = mesh.faces_begin(); it_f != mesh.faces_end(); ++it_f) {
            if (mesh.valence(*it_f) == 3) continue;

            if (mesh.valence(*it_f) == 4) {
                // quad
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
        spdlog::debug("Finished re-mesh to triangle mesh vertices={} faces={} ({:.3}s)",
                      mesh.n_vertices(), mesh.n_faces(), watch);

        return mesh;
    }

    /**
     * Re-meshes a given triangle mesh to a regular isotropic quad mesh.
     * @param sdfn signed distance function to guide the re-meshing
     * @param mesh that is re-meshed
     * @param face_count target face count
     * @return The re-meshed quad mesh.
     */
    entities::Mesh remesh_to_quadmesh(
        const entities::SDFn &sdfn,
        const entities::Mesh &mesh,
        const int face_count
    ) {
        assert(is_trimesh(mesh));
        spdlog::stopwatch watch;

        quadriflow::Parametrizer field;
        spdlog::info("Re-meshing mesh with {} target faces", face_count);


        spdlog::debug("Re-meshing of vertices={}, faces={}", mesh.n_vertices(), mesh.n_faces());

        initialize_parameterizer(field, mesh);
        field.initialize(sdfn, face_count, false);

        spdlog::debug("Finished initializing parameters ({:.3}s)", watch);

        watch.reset();
        spdlog::info("Solving orientation field");

        quadriflow::Optimizer::optimize_orientations(field.m_hierarchy);
        find_orientation_singularities(field.m_hierarchy);

        spdlog::debug("Finished solving orientation field ({:.3}s)", watch);

        watch.reset();
        spdlog::info("Solving field for adaptive scale");

        quadriflow::Optimizer::optimize_scale(field.m_hierarchy, field.m_rho, false);

        spdlog::debug("Finished solving field for adaptive scale ({:.3}s)", watch);


        watch.reset();
        spdlog::info("Solving for position field");

        quadriflow::Optimizer::optimize_positions(field.m_hierarchy);
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

        return from_parametrizer_to_quad_mesh(field);
    }
}
