#pragma once

#include "persistence.h"
#include "services.h"
#include "spdlog/spdlog.h"

namespace bootstrap {
    class Container {
    public:
        Container() {
#ifdef DEV_DEBUG
#include <filesystem>

            spdlog::set_level(spdlog::level::debug);
            if (!std::filesystem::exists("../tests/out/benchmark")) {
                std::filesystem::create_directories("../tests/out/benchmark");
            }
            if (!std::filesystem::exists("../tests/out/stage")) {
                std::filesystem::create_directories("../tests/out/stage");
            }
            if (!std::filesystem::exists("../tests/out/e2e")) {
                std::filesystem::create_directories("../tests/out/e2e");
            }
#endif
        }

        persistence::MeshDao &mesh_dao();

        services::MeshService &mesh_service();

    private:
        std::unique_ptr<persistence::MeshDao> m_mesh_dao = nullptr;
        std::unique_ptr<services::MeshService> m_mesh_service = nullptr;
    };
}
