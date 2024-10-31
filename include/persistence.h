#pragma once

#include <string>
#include "entities.h"

namespace persistence {
    class MeshDao {
    public:
        entities::QuadMesh load_mesh_from_file(const std::string& filename) const;
        void save_mesh_to_file(const std::string& filename, const entities::QuadMesh& mesh) const;
    };
}