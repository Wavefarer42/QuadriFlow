#pragma once

#include "persistence.h"
#include "services.h"
#include "spdlog/spdlog.h"

namespace bootstrap {
    class Container {
    public:
        Container() {
#ifdef DEV_DEBUG
            spdlog::set_level(spdlog::level::debug);
#endif
        }

        persistence::MeshDao &mesh_dao();

        services::MeshService &mesh_service();

    private:
        std::unique_ptr<persistence::MeshDao> m_mesh_dao = nullptr;
        std::unique_ptr<services::MeshService> m_mesh_service = nullptr;
    };
}