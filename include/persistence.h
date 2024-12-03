#pragma once

#include "entities.h"

namespace persistence {
    using namespace Eigen;

    class MeshDao {
    public:
        entities::Mesh load_mesh(const std::string &filename) const;

        void save_mesh(const std::string &filename, const entities::Mesh &mesh) const;

        entities::UnboundModel load_model(const std::string &filename) const;
    };
}