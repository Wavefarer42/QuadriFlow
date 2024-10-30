#pragma once

#include "entities.h"
#include "services.h"

namespace adapters {
    void initialize_parameterizer(services::Parametrizer &field, entities::QuadMesh mesh);
    entities::QuadMesh from_parametrizer_to_quad_mesh(services::Parametrizer &field);
}