#pragma once

#include "entities.h"
#include "services.h"

namespace adapters {

    entities::QuadMesh from_parametrizer_to_quad_mesh(services::Parametrizer &field);
}