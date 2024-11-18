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
}