#pragma once

#include <string>
#include "entities.h"

namespace persistence {
    class MeshDao {
    public:
        static entities::QuadMesh load_mesh_from_file(const std::string& filename);
    };
}