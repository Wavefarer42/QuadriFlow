#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>

#include "gtest/gtest.h"
#include "bootstrap.h"
#include "meshing.h"
#include "sdfn.h"

TEST(LibraryTests, QuadMeshLoading) {
    const persistence::MeshDao sut = bootstrap::Container().mesh_dao();

    const auto mesh = sut.load_mesh("../tests/resources/fandisk.obj");

    for (auto it_f = mesh.faces_begin(); it_f != mesh.faces_end(); ++it_f) {
        EXPECT_EQ(3, mesh.valence(*it_f));
    }
}

TEST(LibraryTests, EigenBroadcasting) {
    Eigen::MatrixXf lhs(4, 3); // Example: 4 rows, 3 columns
    lhs << 1, 2, 3,
            4, 5, 6,
            7, 8, 9,
            10, 11, 12;

    Eigen::VectorXf rhs(4); // Example vector with 4 elements (one per row)
    rhs << 1, 2, 3, 4;

    // Divide each row of the matrix by the corresponding entry in the vector
    Eigen::MatrixXf result = lhs.array().colwise() / rhs.array();
    Eigen::MatrixXf expected = Eigen::MatrixXf(4, 3);
    expected << 1, 2, 3,
            2, 2.5, 3,
            2.33, 2.66, 3,
            2.5, 2.75, 3;

    std::cout << "Result:\n" << result << std::endl;

    EXPECT_TRUE(result.isApprox(expected, 1e-2));
}

TEST(LibraryTests, EigenBroadcastingMultiplication) {
    Eigen::MatrixXf lhs(4, 3); // Example: 4 rows, 3 columns
    lhs << 1, 2, 3,
            1, 2, 3,
            1, 2, 3,
            1, 2, 3;

    Eigen::VectorXf rhs(4); // Example vector with 4 elements (one per row)
    rhs << 1, 2, 3, 4;

    // Divide each row of the matrix by the corresponding entry in the vector
    Eigen::MatrixXf result = lhs.array().colwise() * rhs.array();
    Eigen::MatrixXf expected = Eigen::MatrixXf(4, 3);
    expected << 1, 2, 3,
            2, 4, 6,
            3, 6, 9,
            4, 8, 12;

    std::cout << "Result:\n" << result << std::endl;

    EXPECT_TRUE(result.isApprox(expected, 1e-2));
}

typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;

typedef FT (*Function)(Point_3);

typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

typedef CGAL::Surface_mesh<Point_3> Surface_mesh;


FT sphere_function(Point_3 p) {
    const FT x2 = p.x() * p.x(), y2 = p.y() * p.y(), z2 = p.z() * p.z();
    return x2 + y2 + z2 - 1;
}


TEST(LibraryTests, CGALMeshing) {
    auto model = bootstrap::Container()
            .mesh_dao()
            .load_model("../tests/resources/box-sharp-round-rotated.ubs");

    const auto _sdfn = model[0];
    const auto sdfn = [_sdfn](const Point_3 &p) {
        const Eigen::Vector3f domain = {
            static_cast<float>(p.x()),
            static_cast<float>(p.y()),
            static_cast<float>(p.z())
        };
        const auto distance = _sdfn(domain.transpose());
        return distance(0);
    };

    const auto bounding_box = model.bounding_box(0);
    const auto radius = (bounding_box.max() - bounding_box.min()).norm() / 2;

    Tr tr;
    C2t3 c2t3(tr);
    Surface_3 surface(sdfn, Sphere_3(CGAL::ORIGIN, radius * radius));

    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30., radius / 50, radius / 100);
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
    Surface_mesh sm;
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);
    std::ofstream out("../tests/out/sphere.off");
    out << sm << std::endl;
}
