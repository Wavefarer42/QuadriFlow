#include "spdlog/spdlog.h"

#include "adapters.h"

namespace adapters {

    entities::QuadMesh from_parametrizer_to_quad_mesh(services::Parametrizer &field) {
        spdlog::info("Converting parametrizer to quad mesh");

        entities::QuadMesh mesh_out;

        std::vector<entities::QuadMesh::VertexHandle> handles(field.O_compact.size());
        for (int i = 0; i < field.O_compact.size(); ++i) {
            auto t = field.O_compact[i] * field.normalize_scale + field.normalize_offset;
            handles.emplace_back(
                    mesh_out.add_vertex(entities::QuadMesh::Point(t[0], t[1], t[2]))
            );
        }

        const auto faces_compact = field.F_compact;
        for (const auto &i: faces_compact) {
            mesh_out.add_face({
                                      handles[i[0]],
                                      handles[i[1]],
                                      handles[i[2]],
                                      handles[i[3]]
                              });
        }

        return mesh_out;
    }
}