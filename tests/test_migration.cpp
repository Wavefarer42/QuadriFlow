#include "gtest/gtest.h"

#include "bootstrap.h"
#include "mathext.h"
#include "services.h"

using namespace Eigen;

void normalize_mesh_old(
    MatrixXd &m_vertices,
    double &m_normalize_scale,
    Vector3d &m_normalize_offset
) {
    double maxV[3] = {-1e30, -1e30, -1e30};
    double minV[3] = {1e30, 1e30, 1e30};

    for (int i = 0; i < m_vertices.cols(); ++i) {
        for (int j = 0; j < 3; ++j) {
            maxV[j] = std::max(maxV[j], m_vertices(j, i));
            minV[j] = std::min(minV[j], m_vertices(j, i));
        }
    }
    double scale = std::max(std::max(maxV[0] - minV[0], maxV[1] - minV[1]), maxV[2] - minV[2]) * 0.5;
    for (int i = 0; i < m_vertices.cols(); ++i) {
        for (int j = 0; j < 3; ++j) {
            m_vertices(j, i) = (m_vertices(j, i) - (maxV[j] + minV[j]) * 0.5) / scale;
        }
    }
    m_normalize_scale = scale;
    m_normalize_offset = Vector3d(0.5 * (maxV[0] + minV[0]),
                                  0.5 * (maxV[1] + minV[1]),
                                  0.5 * (maxV[2] + minV[2]));
}


TEST(MigrationSuite, NormalizeMesh) {
    const auto mesh = bootstrap::Container().mesh_service().load_mesh("../tests/resources/box.ply");

    MatrixXd m_vertices(3, mesh.n_vertices());
    for (int i = 0; i < mesh.n_vertices(); ++i) {
        auto v = mesh.point(entities::Mesh::VertexHandle(i));
        m_vertices.col(i) = Vector3d(v[0], v[1], v[2]);
    }

    std::cout << "Original:\n" << m_vertices.block(0, 0, 3, 10) << std::endl;

    double m_normalize_scale = 0;
    Vector3d m_normalize_offset = Vector3d::Zero();
    normalize_mesh_old(m_vertices, m_normalize_scale, m_normalize_offset);

    std::cout << "Normalized:\n" << m_vertices.block(0, 0, 3, 10) << std::endl;

    EXPECT_NE(m_normalize_scale, 0);
    EXPECT_NE(m_normalize_offset, Vector3d::Zero());
}

TEST(MigrationSuite, NormalizeMeshComparison) {
    const auto mesh = bootstrap::Container().mesh_service().load_mesh("../tests/resources/box.ply");

    MatrixXd vertices_orig(3, mesh.n_vertices());
    for (int i = 0; i < mesh.n_vertices(); ++i) {
        auto v = mesh.point(entities::Mesh::VertexHandle(i));
        vertices_orig.col(i) = Vector3d(v[0], v[1], v[2]);
    }

    MatrixXd m_vertices = vertices_orig;
    double m_normalize_scale = 0;
    Vector3d m_normalize_offset = Vector3d::Zero();
    normalize_mesh_old(m_vertices, m_normalize_scale, m_normalize_offset);

    const MatrixXf vertices = vertices_orig.cast<float>().transpose();
    const auto [vertices_scaled, scale, offset] = mathext::normalize(vertices);

    std::cout << "Original:\n" << vertices_orig.block(0, 0, 3, 10) << std::endl;
    std::cout << "Normalized (old):\n" << m_vertices.block(0, 0, 3, 10) << std::endl;
    std::cout << "Normalized (new):\n" << vertices_scaled.block(0, 0, 10, 3).transpose() << std::endl;

    EXPECT_NEAR(m_normalize_scale, scale, 1e-4);
    EXPECT_TRUE(m_normalize_offset.isApprox(offset.cast<double>(), 1e-4));
    EXPECT_TRUE(m_vertices.isApprox(vertices_scaled.cast<double>().transpose(), 1e-4));
}
