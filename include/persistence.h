#pragma once

#include <string>
#include "entities.h"

namespace persistence {
    using namespace Eigen;

    class MeshDao {
    public:
        entities::QuadMesh load_mesh_from_file(const std::string &filename) const;

        void save_mesh_to_file(const std::string &filename, const entities::QuadMesh &mesh) const;

        [[nodiscard]] std::vector<entities::SDFn> load_sdfn_from_file(const std::string &path_model) const;
    };
}