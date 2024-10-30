#pragma once

#include "persistence.h"
#include "services.h"

namespace bootstrap {
    class Container {
    public:
        persistence::MeshDao &mesh_dao();

        services::MeshService &mesh_service();

    private:
        std::unique_ptr<persistence::MeshDao> m_mesh_dao = nullptr;
        std::unique_ptr<services::MeshService> m_mesh_service = nullptr;
    };
}