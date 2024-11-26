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

    entities::Mesh from_parametrizer_to_quad_mesh(const services::Parametrizer &field) {
        spdlog::info("Converting parametrizer to mesh");

        entities::Mesh mesh;

        for (int i = 0; i < field.m_positions_compact.size(); ++i) {
            const Vector3f t = (field.m_positions_compact[i] * field.m_normalize_scale + field.m_normalize_offset).cast<
                float>();
            mesh.add_vertex(entities::Mesh::Point(t[0], t[1], t[2]));
        }

        for (const auto &i: field.m_faces_compact) {
            mesh.add_face({
                entities::Mesh::VertexHandle(i[0]),
                entities::Mesh::VertexHandle(i[1]),
                entities::Mesh::VertexHandle(i[2]),
                entities::Mesh::VertexHandle(i[3])
            });
        }

        return mesh;
    }
}
