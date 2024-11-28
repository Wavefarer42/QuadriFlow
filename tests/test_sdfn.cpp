#include "gtest/gtest.h"

#include "sdfn.h"
#include <Eigen/Core>

#include "bootstrap.h"
#include "mathext.h"

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

TEST(UnboundSuite, BoundingBox) {
    const auto service = bootstrap::Container().mesh_service();
    auto model = service.load_unbound_model_from_file("../tests/resources/box-sharp-aligned.ubs");

    const AlignedBox3f result = model.bounding_box(0);

    EXPECT_EQ(result.min(), Vector3f(-50, -49, -50));
    EXPECT_EQ(result.max(), Vector3f(50, 51, 50));
}

TEST(UnboundSuite, SDFnScale) {
    const auto service = bootstrap::Container().mesh_service();
    auto model = service.load_unbound_model_from_file("../tests/resources/box-sharp-aligned.ubs");


    MatrixXf domain = MatrixXf(4, 3);
    domain << 0, 0, 0,
            0, 1, 0,
            1, 0, 0,
            1, 1, 0;
    const auto [normalized, scale, offset] = mathext::normalize(domain);
    auto sdfn_domain = model[0];
    auto sdfn_normalized = sdfn::scale(sdfn_domain, scale, offset);

    const auto result_domain = sdfn_domain(domain);
    const auto result_normalized = sdfn_normalized(normalized);

    std::cout << "Domain:\n" << result_domain << std::endl;
    std::cout << "Normalized:\n" << result_normalized << std::endl;
    std::cout << "Domain (normalized):\n" << result_normalized * scale << std::endl;

    EXPECT_TRUE(result_domain.isApprox(result_normalized * scale, 1e-4));
}
