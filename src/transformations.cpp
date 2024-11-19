#include <Eigen/Dense>

#include "transformations.h"

namespace transformations {
    Vector4f intersect_planes(
            const MatrixXf &vertices,
            const MatrixXf &normals
    ) {
        const auto svd_vmul_sym = [](MatrixXf a, VectorXf v) {
            return Vector3f{
                    a.row(0).dot(v),
                    a(0, 1) * v[0] + a(1, 1) * v[1] + a(1, 2) * v[2],
                    a(0, 2) * v[0] + a(1, 2) * v[1] + a(2, 2) * v[2]
            };
        };

        const MatrixXf ATA = normals.transpose() * normals;
        const VectorXf b = (vertices.array() * normals.array()).rowwise().sum();
        VectorXf ATb = (b.asDiagonal() * normals).colwise().sum();
        const Vector3f centroid = normals.colwise().mean() / static_cast<float>(vertices.rows());

        const VectorXf ATb_mass = svd_vmul_sym(ATA, centroid);
        ATb -= ATb_mass.transpose();

        VectorXf x = ATA.fullPivLu().solve(ATb);
        const Vector3f differences = ATb - svd_vmul_sym(ATA, x);
        float error = differences.dot(differences);

        x += centroid;

        return {x[0], x[1], x[2], error};
    }
}