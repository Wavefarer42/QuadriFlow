#pragma once

#include <string>
#include "entities.h"

namespace persistence {
    using namespace Eigen;

    class MeshDao {
    public:
        entities::Mesh load_mesh_from_file(const std::string &filename) const;

        void save_mesh_to_file(const std::string &filename, const entities::Mesh &mesh) const;

        entities::UnboundModel load_unbound_model(const std::string &filename) const;
    };
}