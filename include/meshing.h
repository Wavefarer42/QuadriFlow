#pragma once

#include "entities.h"

namespace surfacenets {
    using namespace Eigen;
    using NdToFlatIndexer = std::function<int(const Vector3i &)>;

    NdToFlatIndexer indexer_nd_to_linear(
        int resolution
    );

    MatrixXi create_indices(
        int resolution
    );

    MatrixXf scale_to_domain(
        const MatrixXi &indices,
        const AlignedBox3f &bounds,
        const int resolution
    );

    MatrixXf sample_sdf(
        const entities::SDFn &sdfn,
        const MatrixXf &domain
    );

    std::vector<Vector3f> create_vertices(
        MatrixXi &indices,
        const MatrixXf &domain,
        const MatrixXf &sdf,
        const int resolution,
        const NdToFlatIndexer &linearize
    );
}

namespace meshing {
    using namespace Eigen;

    bool is_trimesh(
        const entities::Mesh &mesh
    );

    entities::Mesh mesh_to_trimesh(
        const entities::SDFn &sdfn,
        const AlignedBox3f &bounds,
        int resolution = 100,
        const std::string &algorithm = "delaunay"
    );

    entities::Mesh mesh_to_quadmesh(
        const entities::SDFn &sdfn,
        const AlignedBox3f &bounds,
        int resolution = 100,
        const std::string &algorithm = "surface-nets"
    );

    entities::Mesh remesh_to_trimesh(
        entities::Mesh &mesh
    );

    entities::Mesh remesh_to_quadmesh(
        const entities::SDFn &sdfn,
        const entities::Mesh &mesh,
        int face_count
    );
};
