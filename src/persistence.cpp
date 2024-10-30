#include "persistence.h"
#include <OpenMesh/Core/IO/MeshIO.hh>

namespace persistence {
    entities::QuadMesh MeshDao::load_mesh_from_file(const std::string &filename) {

        entities::QuadMesh mesh;
        if (!OpenMesh::IO::read_mesh(mesh, filename)) {
            throw std::runtime_error("Could not read mesh from file" + filename);
        }

        return mesh;
    }
}