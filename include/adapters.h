#pragma once

#include "entities.h"
#include "services.h"

namespace adapters {
    void initialize_parameterizer(services::Parametrizer &field, entities::Mesh mesh);
    entities::Mesh from_parametrizer_to_quad_mesh(services::Parametrizer &field);
}