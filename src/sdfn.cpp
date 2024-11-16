#include "sdfn.h"

namespace sdfn {
    VectorXf sphere(const MatrixXf &domain) {
        float radius = 1;
        const Vector3f origin = Vector3f::Zero();
        return (domain.rowwise() - origin.transpose()).rowwise().norm().array() - radius;
    }

    VectorXf box(const MatrixXf domain) {
        const Vector3f half_size(0.5f, 0.5f, 0.5f);
        VectorXf distances(domain.rows());
        tbb::parallel_for(
                tbb::blocked_range<int>(0, domain.rows()),
                [&](const tbb::blocked_range<int> &range) {
                    for (int i = range.begin(); i < range.end(); ++i) {
                        Vector3f point = domain.row(i);

                        Vector3f d = (point.cwiseAbs() - half_size).cwiseMax(0.0f);
                        float outside_distance = d.norm();
                        float inside_distance = std::min(std::max(point.x(), -half_size.x()), half_size.x());
                        inside_distance = std::min(inside_distance, std::max(point.y(), -half_size.y()));
                        inside_distance = std::min(inside_distance, std::max(point.z(), -half_size.z()));

                        distances(i) = (point.cwiseAbs().maxCoeff() <= half_size.maxCoeff())
                                       ? inside_distance : outside_distance;
                    }
                }
        );

        return distances;
    }
}