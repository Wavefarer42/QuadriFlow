#pragma once

#include "entities.h"

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
