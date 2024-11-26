#include "spdlog/spdlog.h"

#include "adapters.h"

namespace adapters {
    using namespace Eigen;

    void initialize_parameterizer(services::Parametrizer &field, entities::Mesh mesh) {
        spdlog::debug("Initializing parameters");

        field.m_vertices = MatrixXd(3, mesh.n_vertices());
        for (auto it_v = mesh.vertices_begin(); it_v != mesh.vertices_end(); ++it_v) {
            auto idx = (*it_v).idx();
            auto point = mesh.point(*it_v);
            field.m_vertices.col(idx) = Vector3d(point[0], point[1], point[2]);
        }

        field.m_faces = MatrixXi(3, mesh.n_faces());
        for (auto it_f = mesh.faces_begin(); it_f != mesh.faces_end(); ++it_f) {
            auto idx = (*it_f).idx();
            auto fv_it = mesh.cfv_iter(*it_f);
            for (int i = 0; i < 3; ++i) {
                field.m_faces(i, idx) = (*fv_it).idx();
                ++fv_it;
            }
        }
    }

    entities::Mesh from_parametrizer_to_quad_mesh(services::Parametrizer &field) {
        spdlog::info("Converting parametrizer to quad mesh");

        entities::Mesh mesh_out;

        std::vector<entities::Mesh::VertexHandle> handles(field.m_positions_compact.size());
        for (int i = 0; i < field.m_positions_compact.size(); ++i) {
            auto t = field.m_positions_compact[i] * field.m_normalize_scale + field.m_normalize_offset;
            handles.emplace_back(
                    mesh_out.add_vertex(entities::Mesh::Point(t[0], t[1], t[2]))
            );
        }

        const auto faces_compact = field.m_faces_compact;
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