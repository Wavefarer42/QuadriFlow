#include <fstream>
#include <stdexcept>
#include <boost/process.hpp>
#include <nlohmann/json.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "spdlog/spdlog.h"

#include "persistence.h"


namespace persistence {
    MatrixXf from_list_list_to_matrix(std::vector<std::vector<float> > data) {
        auto rows = data.size();
        auto cols = data[0].size();
        MatrixXf mat(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                mat(i, j) = data[i][j];
            }
        }
        return mat;
    }

    MatrixXf from_list_to_vector(std::vector<float> data) {
        auto rows = data.size();
        MatrixXf mat(rows, 1);
        for (int i = 0; i < rows; ++i) {
            mat(i, 0) = data[i];
        }
        return mat;
    }

    std::vector<std::vector<float> > from_matrix_to_list_list(MatrixXf data) {
        auto rows = data.rows();
        auto cols = data.cols();
        std::vector<std::vector<float> > vec(rows, std::vector<float>(cols));
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                vec[i][j] = data(i, j);
            }
        }
        return vec;
    }

    entities::Mesh MeshDao::load_mesh(const std::string &filename) const {
        entities::Mesh mesh;
        if (!OpenMesh::IO::read_mesh(mesh, filename)) {
            throw std::runtime_error("Could not read mesh from file " + filename);
        }

        return mesh;
    }

    void MeshDao::save_mesh(
        const std::string &filename,
        const entities::Mesh &mesh
    ) const {
        if (!OpenMesh::IO::write_mesh(mesh, filename)) {
            throw std::runtime_error("Could not write mesh to file " + filename);
        }
    }

    entities::UnboundModel MeshDao::load_model(const std::string &path_model) const {
        assert(path_model.ends_with(".ubs"));
        spdlog::info("Loading unbound model for {}", path_model);


        std::ifstream f(path_model);
        if (f.is_open()) {
            try {
                nlohmann::json data = nlohmann::json::parse(f);
                const auto collection = data.get<UB::Collection>();
                return entities::UnboundModel(collection);
            } catch (nlohmann::json::exception const &e) {
                spdlog::error("Could not parse collection json:\n{}", e.what());
                throw std::ios_base::failure("Could not parse collection json");
            }
        } else {
            spdlog::error("Could not open the collection file.");
            throw std::ios_base::failure("Could not open the collection file.");
        }
    }
}
