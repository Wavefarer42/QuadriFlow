#include "gtest/gtest.h"

#include "sdfn.h"
#include "bootstrap.h"
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace Eigen;

TEST(E2E, GradientOfSphere) {
    auto domain = MatrixXf(4, 3);
    domain << 0, 0, 0,
            0, 1, 0,
            1, 0, 0,
            1, 1, 0;

    const MatrixXf result = sdfn::gradient_of(sdfn::sphere, domain);

    MatrixXf expected = MatrixXf(domain.rows(), domain.cols());
    expected << 0, 0, 0,
            0, 1, 0,
            1, 0, 0,
            0.70, 0.70, 0;

    std::cout << result << std::endl;
    ASSERT_TRUE(result.isApprox(expected, 1e-2));
}
