#pragma once

#include "entities.h"
#include "persistence.h"

namespace services {
    class MeshService {
    public:
        explicit MeshService(persistence::MeshDao mesh_dao) : mesh_dao(mesh_dao) {
        }

        void to_isotropic_quadmesh(
            const std::string &path_input,
            const std::string &path_output,
            int face_count = 10000,
            int sdfn_resolution = 100
        ) const;

    private:
        persistence::MeshDao mesh_dao;
    };
}
