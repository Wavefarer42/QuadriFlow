#include <fstream>

#include <boost/process.hpp>
#include <nlohmann/json.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "spdlog/spdlog.h"
#include "persistence.h"

#ifndef PATH_UBTOOL
#define PATH_UBTOOL "UbTool"
#endif


namespace persistence {

    MatrixXf from_list_list_to_matrix(std::vector<std::vector<float>> data) {
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

    std::vector<std::vector<float>> from_matrix_to_list_list(MatrixXf data) {
        auto rows = data.rows();
        auto cols = data.cols();
        std::vector<std::vector<float>> vec(rows, std::vector<float>(cols));
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                vec[i][j] = data(i, j);
            }
        }
        return vec;
    }

    entities::QuadMesh MeshDao::load_mesh_from_file(const std::string &filename) const {

        entities::QuadMesh mesh;
        if (!OpenMesh::IO::read_mesh(mesh, filename)) {
            throw std::runtime_error("Could not read mesh from file " + filename);
        }

        return mesh;
    }

    void MeshDao::save_mesh_to_file(const std::string &filename,
                                    const entities::QuadMesh &mesh) const {
        if (!OpenMesh::IO::write_mesh(mesh, filename)) {
            throw std::runtime_error("Could not write mesh to file " + filename);
        }
    }

    std::vector<entities::UnboundModelInfo> load_model_info(const std::string &path_model) {
        spdlog::info("Loading model info for {}", path_model);

        const auto tool = std::string(PATH_UBTOOL);
        const auto command = std::format("{} info -u {}", tool, path_model);

        boost::process::ipstream pipe_stream;
        boost::process::child process(command, boost::process::std_out > pipe_stream);

        std::string json_data;

        std::string line;
        while (pipe_stream && std::getline(pipe_stream, line) && !line.empty()) {
            json_data += line;
        }

        process.wait();

        if (process.exit_code() != 0 || json_data.empty()) {
            throw std::runtime_error("Expected the model infos for " + path_model);
        }

        nlohmann::json json = nlohmann::json::parse(json_data);

        std::vector<entities::UnboundModelInfo> infos;
        for (auto &[id_model, entry]: json.items()) {
            infos.push_back(
                    {
                            .id_model = id_model,
                            .size_edits = entry["size_edits"],
                            .resolution = entry["resolution"],
                            .shape_density = entry["shape_density"]
                    }
            );
        }

        return infos;
    }

    void write_sample_domain(
            const std::string &path_data,
            const std::string &id_model,
            const MatrixXf &domain
    ) {
        spdlog::debug("Writing sample domain data to {}", path_data);

        nlohmann::json json;
        json[id_model] = from_matrix_to_list_list(domain);

        std::ofstream f(path_data);
        if (f.is_open()) {
            f << json.dump(4);
            f.close();
        } else {
            spdlog::error("Could not open file {}", path_data);
            throw std::runtime_error("Could not write sample domain data " + path_data);
        }
    }

    MatrixXf sample_model(
            const std::string &path_model,
            const std::string &id_model,
            const MatrixXf &domain
    ) {
        spdlog::debug("Loading model sample for {}", path_model);

        const auto path_domain = "domain.json";
        write_sample_domain(path_domain, id_model, domain);

        const auto tool = std::string(PATH_UBTOOL);
        const auto command = std::format("{} sample -u {} -d {}", tool, path_model, path_domain);

        boost::process::ipstream pipe_stream;
        boost::process::child process(command, boost::process::std_out > pipe_stream);

        std::string json_data;

        std::string line;
        while (pipe_stream && std::getline(pipe_stream, line) && !line.empty()) {
            json_data += line;
        }

        process.wait();

        if (process.exit_code() != 0 || json_data.empty()) {
            throw std::runtime_error("Expected the model infos for " + path_model);
        }

        nlohmann::json json = nlohmann::json::parse(json_data);
        if (!json.contains(id_model)) {
            throw std::runtime_error("Expected the model sample in the UbTool response for " + id_model);
        }

        const auto data = json[id_model];
        const auto distances = from_list_to_vector(data["distances"]);
        const auto gradients = from_list_list_to_matrix(data["gradients"]);

        MatrixXf result(domain.rows(), 4);
        result << distances, gradients;

        return result;
    }

    std::vector<entities::SDFn> MeshDao::load_sdfn_from_file(const std::string &path_model) const {
        spdlog::debug("Creating new Unbound sampler for file {}", path_model);

        const auto models = load_model_info(path_model);

        std::vector<entities::SDFn> sdfns;
        for (const auto &model: models) {
            const auto sampler = [model, path_model](const MatrixXf &domain) {
                return sample_model(path_model, model.id_model, domain);
            };

            sdfns.emplace_back(sampler);
        }

        return sdfns;
    }
}