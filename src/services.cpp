#include <OpenMesh/Core/IO/MeshIO.hh>
#include "services.h"

namespace services {

    entities::QuadMesh MeshService::load_trimesh_from_file(const std::string &filename) {
        const auto mesh = this->mesh_dao.load_mesh_from_file(filename);

        // Validate that the mesh is a triangle mesh
        int n_triangles = 0;
        int n_non_triangles = 0;
        for (auto it_f = mesh.faces_begin(); it_f != mesh.faces_end(); ++it_f) {
            if (mesh.valence(*it_f) == 3) {
                n_triangles++;
            } else {
                n_non_triangles++;
            }
        }

        if (n_non_triangles > 0) {
            throw std::runtime_error(
                    std::format("Please provide a triangle mesh as input. Triangles={}, Non Triangle={}",
                                n_triangles, n_non_triangles)
            );
        }

        return mesh;
    }

    void MeshService::save_quadmesh_to_file(const std::string &filename, const entities::QuadMesh &mesh) {
        this->mesh_dao.save_mesh_to_file(filename, mesh);
    }
}