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
    static const std::vector AXES = {AXIS_X, AXIS_Y, AXIS_Z};

    static const MatrixXi QUAD_POINTS_X = (Matrix<int, 4, 3>() <<
                                           0, 0, -1,
                                           0, -1, -1,
                                           0, -1, 0,
                                           0, 0, 0).finished();
    static const MatrixXi QUAD_POINTS_Y = (Matrix<int, 4, 3>() <<
                                           0, 0, -1,
                                           0, 0, 0,
                                           -1, 0, 0,
                                           -1, 0, -1).finished();
    static const MatrixXi QUAD_POINTS_Z = (Matrix<int, 4, 3>() <<
                                           0, 0, 0,
                                           0, -1, 0,
                                           -1, -1, 0,
                                           -1, 0, 0).finished();
    static const std::vector QUAD_POINTS = {QUAD_POINTS_X, QUAD_POINTS_Y, QUAD_POINTS_Z};

    Vector3f estimate_centroid(
        const VectorXf &distances_corners
    ) {
        int n = 0;
        Vector3f centroid = Vector3f::Zero();
        for (int i = 0; i < CUBE_EDGES.rows(); ++i) {
            const auto idx_edge_a = CUBE_EDGES(i, 0);
            const auto idx_edge_b = CUBE_EDGES(i, 1);
            const auto distance_a = distances_corners(idx_edge_a);
            const auto distance_b = distances_corners(idx_edge_b);

            if (distance_a < 0.0f != distance_b < 0.0f) {
                const auto interpolation1 = distance_a / (distance_a - distance_b + 1e-10f);
                const auto interpolation2 = 1.0f - interpolation1;

                const Vector3f interpolation = interpolation2 * CUBE_CORNERS.row(idx_edge_a).cast<float>()
                                               + interpolation1 * CUBE_CORNERS.row(idx_edge_b).cast<float>();
                centroid += interpolation;
                n++;
            }
        }

        assert(n > 0);
        return centroid / static_cast<float>(n);
    }

    MatrixXi create_indices(
        int resolution
    ) {
        spdlog::debug("Creating domain indices");
        spdlog::stopwatch watch;

        // [0:3]: Grid coordinate (x, y, z)
        // [3]: Mapping from grid coordinate to linear surface array index see vertices
        const auto grid_max = resolution + 1;
        MatrixXi indices(grid_max * grid_max * grid_max, 3);
        int index = 0;
        for (int z = 0; z < grid_max; ++z) {
            for (int y = 0; y < grid_max; ++y) {
                for (int x = 0; x < grid_max; ++x) {
                    indices.row(index++) << x, y, z;
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

    bool is_on_surface(
        const float distance1,
        const float distance2
    ) {
        return (distance1 < 0.0f && distance2 >= 0.0f) || (distance2 < 0.0f && distance1 >= 0.0f);
    }

    bool is_negative_face(
        const float distance1,
        const float distance2
    ) {
        return distance2 < 0.0f && distance1 >= 0.0f;
    }

    entities::Mesh create_mesh(
        MatrixXi &indices,
        const MatrixXf &domain,
        const MatrixXf &sdf,
        const int resolution,
        const NdToFlatIndexer &linearize
    ) {
        spdlog::debug("Creating surface mesh");
        spdlog::stopwatch watch;

        entities::Mesh mesh;

        const Vector3f domain_offset = domain.row(0);
        const Vector3f domain_scale = (domain.row(domain.rows() - 1) - domain.row(0)).array() / resolution;
        std::unordered_map<int, entities::Mesh::VertexHandle> domain_to_surface;

        for (int idx_cell = 0; idx_cell < indices.rows(); ++idx_cell) {
            if (indices(idx_cell, 0) >= resolution
                || indices(idx_cell, 1) >= resolution
                || indices(idx_cell, 2) >= resolution) {
                continue;
            }
            const Vector3i corner = indices.row(idx_cell).head<3>();

            const std::set<int> debugstop = {
                31789, 31739, 31790,
                29189, 31789, 31790,
                26588, 29138, 29189,
                29138, 31688, 31739
            };
            if (debugstop.contains(idx_cell)) {
                spdlog::debug("Vertex at  with distances=");
            }

            for (int idx_axis = 0; idx_axis < AXES.size(); ++idx_axis) {
                const auto &axis = AXES[idx_axis];

                const int idx1 = linearize(corner);
                const int idx2 = linearize(corner + axis);
                const auto d1 = sdf(idx1);
                const auto d2 = sdf(idx2);

                if (is_on_surface(d1, d2)) {
                    std::vector<entities::Mesh::VertexHandle> face(4);

                    // Create vertex if needed
                    VectorXf distances_corners(CUBE_CORNERS.rows());
                    for (int j = 0; j < CUBE_CORNERS.rows(); ++j) {
                        const Vector3i idx_nd = corner + CUBE_CORNERS.row(j).transpose();
                        const int idx_flat = linearize(idx_nd);
                        distances_corners(j) = sdf(idx_flat);
                    }

                    std::cout << "corner=" << corner.transpose() << std::endl << " distances=" << distances_corners.
                            transpose()
                            << std::endl;

                    const MatrixXi indices_quads = QUAD_POINTS[idx_axis].rowwise() + corner.transpose();
                    std::cout << "quads points=" << indices_quads << std::endl;

                    for (int idx_vertex = 0; idx_vertex < indices_quads.rows(); ++idx_vertex) {
                        const auto idx_vertex_flat = linearize(indices_quads.row(idx_vertex));

                        if (!domain_to_surface.contains(idx_vertex_flat)) {
                            // const Vector3f centroid = corner.cast<float>();
                            const Vector3f centroid = estimate_centroid(distances_corners) + corner.cast<float>();
                            const Vector3f position =
                                    domain_offset + (centroid.array() * domain_scale.array()).matrix();

                            const auto idx_v = mesh.add_vertex(
                                entities::Mesh::Point(position.x(), position.y(), position.z())
                            );

                            const std::set<int> debugstop = {
                                610, 605, 611,
                                491, 610, 611,
                                380, 487, 490,
                                487, 602, 605
                            };
                            if (debugstop.contains(idx_v.idx())) {
                                spdlog::debug("Vertex at  with distances=");
                            }
                            domain_to_surface[idx_vertex_flat] = idx_v;
                        }

                        face[idx_vertex] = domain_to_surface[idx_vertex_flat];
                    }

                    if (is_negative_face(d1, d2)) {
                        const auto fv = mesh.add_face(
                            face[0],
                            face[1],
                            face[2],
                            face[3]
                        );
                    } else {
                        const auto fv = mesh.add_face(
                            face[3],
                            face[2],
                            face[1],
                            face[0]
                        );
                    }
                }
            }
        }

        spdlog::debug("Finished creating surface vertices={} ({:.3}s)", domain_to_surface.size(), watch);

#ifdef DEV_DEBUG
        OpenMesh::IO::write_mesh(mesh, "../tests/out/stage/4-mesh.ply");
#endif

        return mesh;
    }

    /**
     * Creates a index linearizer into the original samples grid index space.
     */
    NdToFlatIndexer indexer_nd_to_linear(
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
        const auto mesh = create_mesh(indices, domain, sdf, resolution, linearize);

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
        const auto distance_bound = diameter / 100;

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
