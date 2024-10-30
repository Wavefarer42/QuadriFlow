#include "bootstrap.h"

namespace bootstrap {
    persistence::MeshDao &Container::mesh_dao() {
        if (m_mesh_dao == nullptr) {
            m_mesh_dao = std::make_unique<persistence::MeshDao>();
        }
        return *m_mesh_dao;
    }

    services::MeshService &Container::mesh_service() {
        if (m_mesh_service == nullptr) {
            m_mesh_service = std::make_unique<services::MeshService>(mesh_dao());
        }
        return *m_mesh_service;
    }
}