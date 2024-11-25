#include "sdfn.h"

namespace sdfn {
    VectorXf sphere(const MatrixXf &domain) {
        float radius = 1;
        const Vector3f origin = Vector3f::Zero();
        return (domain.rowwise() - origin.transpose()).rowwise().norm().array() - radius;
    }


    VectorXf box(const MatrixXf &domain) {
        const Vector3f boxSize = Vector3f(1, 1, 1);
        MatrixXf q = domain.array().cwiseAbs().rowwise() - boxSize.transpose().array();

        MatrixXf a = q.cwiseMax(0);
        VectorXf b = q.rowwise().maxCoeff().cwiseMin(0);

        VectorXf distances = a.rowwise().norm() + b;

        return distances;
    }

    VectorXf cylinder(const MatrixXf &domain) {
        const float radius = 0.8;
        const VectorXf distances = domain.block(0, 1, domain.rows(), 2).rowwise().norm().array() - radius;

        return distances;
    }

    entities::SDFn rotate(
        entities::SDFn sdfn,
        const Vector3f axis,
        const float angle
    ) {
        const auto _rotate = [sdfn, axis, angle](const MatrixXf &domain) -> VectorXf {
            const AngleAxisf rotation(angle, axis.normalized());
            const Matrix3f rotationMatrix = rotation.toRotationMatrix();

            // Apply the inverse rotation to each point in the domain
            const MatrixXf rotatedDomain = (rotationMatrix.inverse() * domain.transpose()).transpose();

            const VectorXf distances = sdfn(rotatedDomain);

            return distances;
        };

        return _rotate;
    }

    MatrixXf gradient_of(
        const entities::SDFn &sdfn,
        const MatrixXf &domain,
        const float epsilon
    ) {
        auto gradient = MatrixXf(domain.rows(), 3);

        const auto eps = MatrixXf::Identity(3, 3) * epsilon;

        gradient.col(0) = (sdfn(domain.rowwise() + eps.row(0)) - sdfn(domain.rowwise() - eps.row(0))).array()
                          / (2 * epsilon);
        gradient.col(1) = (sdfn(domain.rowwise() + eps.row(1)) - sdfn(domain.rowwise() - eps.row(1))).array()
                          / (2 * epsilon);
        gradient.col(2) = (sdfn(domain.rowwise() + eps.row(2)) - sdfn(domain.rowwise() - eps.row(2))).array()
                          / (2 * epsilon);

        return gradient;
    }

    MatrixXf normal_of(
        const entities::SDFn &sdfn,
        const MatrixXf &domain,
        const float epsilon
    ) {
        const MatrixXf gradient = gradient_of(sdfn, domain, epsilon);
        return normal_of(gradient);
    }

    MatrixXf normal_of(
        const MatrixXf &gradients
    ) {
        VectorXf norms = gradients.rowwise().norm();
        norms = (norms.array() == 0).select(1, norms);


        return gradients.array().colwise() / norms.array();
    }
}
