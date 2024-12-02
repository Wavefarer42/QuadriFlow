#include "spdlog/spdlog.h"

#include "mathext.h"

namespace mathext {
    MatrixXf clip(
        const MatrixXf &mat,
        float minVal,
        float maxVal
    ) {
        return mat.unaryExpr([minVal, maxVal](float val) {
            return std::min(std::max(val, minVal), maxVal);
        });
    }

    float frobenius_norm_off_diagonal(const Eigen::MatrixXf &A) {
        int n = A.rows();
        float sum = 0.0f;

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += A(i, j) * A(i, j);
                }
            }
        }

        return std::sqrt(sum);
    }

    float percentile(const Eigen::VectorXf &vec, float percentile) {
        // Copy the data from the Eigen vector to a standard vector for sorting
        Eigen::VectorXf vec_eval = vec.eval();
        std::vector<float> sorted_vec(vec_eval.data(), vec_eval.data() + vec_eval.size());

        // Sort the data
        std::sort(sorted_vec.begin(), sorted_vec.end());

        // Compute the percentile index
        float percentile_index = percentile * (sorted_vec.size() - 1);
        int lower_index = static_cast<int>(percentile_index);
        int upper_index = lower_index + 1;

        // Get the lower and upper values
        float lower_value = sorted_vec[lower_index];
        float upper_value = (upper_index < sorted_vec.size()) ? sorted_vec[upper_index] : lower_value;

        // Linear interpolation if needed
        float weight = percentile_index - lower_index;
        float percentile_value = lower_value * (1.0f - weight) + upper_value * weight;

        return percentile_value;
    }

    // Mesh math
    Vector3f face_centroid(
        entities::Mesh &mesh,
        const entities::Mesh::FaceHandle &face
    ) {
        int total = 0;
        Vector3f centroid = Vector3f::Zero();
        for (auto it_face_vertex = mesh.fv_iter(face);
             it_face_vertex.is_valid(); ++it_face_vertex) {
            const auto p = mesh.point(*it_face_vertex);
            centroid += Vector3f(p[0], p[1], p[2]);
            total++;
        }

        centroid /= total;

        return centroid;
    }

    MatrixXf face_centroids_ring(
        entities::Mesh &mesh,
        const entities::Mesh::VertexHandle vertex
    ) {
        std::vector<Vector3f> centroids_list;
        for (auto it_face = mesh.vf_iter(vertex);
             it_face.is_valid(); ++it_face) {
            const auto centroid = face_centroid(mesh, *it_face);
            centroids_list.push_back(centroid);
        }

        MatrixXf centroids(centroids_list.size(), 3);
        for (int i = 0; i < centroids_list.size(); ++i) {
            centroids.row(i) = centroids_list[i];
        }

        return centroids;
    }

    MatrixXf count_unique(
        const MatrixXf &mat
    ) {
        // Create an unordered map to store unique elements and their counts
        std::map<float, int> elementCounts;

        // Iterate over all elements of the matrix and populate the map with counts
        for (int i = 0; i < mat.rows(); ++i) {
            for (int j = 0; j < mat.cols(); ++j) {
                elementCounts[mat(i, j)]++;
            }
        }

        // Create a MatrixXf with two columns: one for unique elements, one for their counts
        int uniqueCount = elementCounts.size();
        MatrixXf uniqueElementsWithCounts(uniqueCount, 2);

        // Fill the matrix with unique elements and their counts
        int index = 0;
        for (const auto &[element, count]: elementCounts) {
            uniqueElementsWithCounts(index, 0) = element; // Unique element
            uniqueElementsWithCounts(index, 1) = count; // Count of that element
            index++;
        }

        return uniqueElementsWithCounts;
    }

    Vector4f _intersect_planes(
        const MatrixXf &vertices,
        const MatrixXf &normals,
        const Vector3f &bias
    ) {
        const auto svd_vmul_sym = [](MatrixXf a, VectorXf v) {
            return Vector3f{
                a.row(0).dot(v),
                a(0, 1) * v[0] + a(1, 1) * v[1] + a(1, 2) * v[2],
                a(0, 2) * v[0] + a(1, 2) * v[1] + a(2, 2) * v[2]
            };
        };

        const MatrixXf ATA = (normals.transpose() * normals).triangularView<Upper>();
        const VectorXf b = (vertices.array() * normals.array()).rowwise().sum();
        VectorXf ATb = (b.asDiagonal() * normals).colwise().sum();

        const VectorXf constraint_center = svd_vmul_sym(ATA, bias);
        ATb -= constraint_center.transpose();

        VectorXf x = ATA.colPivHouseholderQr().solve(ATb);
        const Vector3f differences = ATb - svd_vmul_sym(ATA, x);
        float error = differences.dot(differences);

        x += bias;

        return {x[0], x[1], x[2], error};
    }

    Vector4f intersect_planes(
        const MatrixXf &vertices,
        const MatrixXf &normals,
        const Vector3f &bias,
        int patience
    ) {
        assert(vertices.rows() == normals.rows());

        const float max_displacement = (vertices.rowwise() - bias.transpose())
                                       .rowwise()
                                       .norm()
                                       .minCoeff() * 0.9;

        float displacement = 0;
        Vector4f solution;
        int iteration = 0;
        do {
            solution = _intersect_planes(vertices, normals, bias);
            displacement = (solution.head<3>() - bias).norm();

            iteration++;
        } while (displacement > max_displacement && iteration < patience);

        return solution;
    }

    Vector4f intersect_planes(
        const MatrixXf &vertices,
        const MatrixXf &normals,
        int patience
    ) {
        const Vector3f centroid = vertices.colwise().sum() / static_cast<float>(vertices.rows());
        return intersect_planes(vertices, normals, centroid, patience);
    }


    std::tuple<MatrixXf, float, Vector3f> normalize(
        const MatrixXf &vertices
    ) {
        const Vector3f upper = vertices.colwise().maxCoeff();
        const Vector3f lower = vertices.colwise().minCoeff();

        const float scale = (upper - lower).maxCoeff() * 0.5;
        const Vector3f offset = (upper + lower) * 0.5;
        const MatrixXf vertices_scaled = (vertices.rowwise() - offset.transpose()) / scale;

        return std::make_tuple(
            vertices_scaled,
            scale,
            offset
        );
    }

    MatrixXf denormalize(
        const MatrixXf &vertices,
        const float scale,
        const Vector3f &offset
    ) {
        return (vertices * scale).rowwise() + offset.transpose();
    }
}
