#include "svd.h"

namespace svd {
    using namespace Eigen;

    std::tuple<MatrixXf, MatrixXf, Vector3f> create_linear_system(
            const MatrixXf &vertices,
            const MatrixXf &normals
    ) {
        const MatrixXf ATA = normals.transpose() * normals;
        const VectorXf b = (vertices.array() * normals.array()).rowwise().sum();
        const MatrixXf ATb = (b.asDiagonal() * normals).colwise().sum();
        const Vector3f centroid = normals.colwise().mean();
        const Vector4f centroid_weight = {
                centroid[0],
                centroid[1],
                centroid[2],
                static_cast<float>(vertices.rows())
        };

        return std::make_tuple(ATA, ATb, centroid_weight);
    }

    Vector3f svd_vmul_sym(MatrixXf a, VectorXf v) {
        return {
                a.row(0).dot(v),
                a(0, 1) * v[0] + a(1, 1) * v[1] + a(1, 2) * v[2],
                a(0, 2) * v[0] + a(1, 2) * v[1] + a(2, 2) * v[2]
        };
    }

    Vector4f intersect_planes(
            const MatrixXf &vertices,
            const MatrixXf &normals
    ) {
        auto [ATA, ATb, centroid] = create_linear_system(vertices, normals);

        const Vector3f masspoint = centroid.head<3>() / centroid[3];
        ATb -= svd_vmul_sym(ATA, masspoint);

        VectorXf x = ATA.llt().solve(ATb);
        const Vector3f differences = ATb - svd_vmul_sym(ATA, x);
        float error = differences.dot(differences);

        x += masspoint;

        return {x[0], x[1], x[2], error};
    }
}
