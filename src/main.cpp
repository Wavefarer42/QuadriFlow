#include <cstdlib>
#include <chrono>
#include <format>

#include "OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh"
#include "adapters.h"
#include "bootstrap.h"
#include "persistence.h"
#include "services.h"
#include "optimizer.h"
#include "field-math.h"

using namespace services;

Parametrizer field;

unsigned long long inline GetCurrentTime64() {
    using namespace std::chrono;
    return duration_cast<milliseconds>(steady_clock::now().time_since_epoch()).count();
}


int main(int argc, char **argv) {
    setbuf(stdout, NULL);

    int t1, t2;
    std::string input_obj, output_obj;
    int faces = -1;
    for (int i = 0; i < argc; ++i) {
        if (strcmp(argv[i], "-f") == 0) {
            sscanf(argv[i + 1], "%d", &faces);
        } else if (strcmp(argv[i], "-i") == 0) {
            input_obj = argv[i + 1];
        } else if (strcmp(argv[i], "-o") == 0) {
            output_obj = argv[i + 1];
        } else if (strcmp(argv[i], "-sharp") == 0) {
            field.flag_preserve_sharp = 1;
        } else if (strcmp(argv[i], "-boundary") == 0) {
            field.flag_preserve_boundary = 1;
        } else if (strcmp(argv[i], "-adaptive") == 0) {
            field.flag_adaptive_scale = 1;
        } else if (strcmp(argv[i], "-seed") == 0) {
            field.m_hierarchy.rng_seed = atoi(argv[i + 1]);
        }
    }
    printf("%d %s %s\n", faces, input_obj.c_str(), output_obj.c_str());

    bootstrap::Container container = bootstrap::Container();
    MeshService service = container.mesh_service();
    const auto mesh = service.load_trimesh_from_file(input_obj);

    field.m_vertices = MatrixXd(3, mesh.n_vertices());
    for (auto it_v = mesh.vertices_begin(); it_v != mesh.vertices_end(); ++it_v) {
        auto idx = (*it_v).idx();
        auto point = mesh.point(*it_v);
        field.m_vertices.col(idx) = Vector3d(point[0], point[1], point[2]);
    }

    field.m_faces = MatrixXi(3, mesh.n_faces());
    for (auto it_f = mesh.faces_begin(); it_f != mesh.faces_end(); ++it_f) {
        auto idx = (*it_f).idx();
        auto fv_it = mesh.cfv_iter(*it_f);
        for (int i = 0; i < 3; ++i) {
            field.m_faces(i, idx) = (*fv_it).idx();
            ++fv_it;
        }
    }

    printf("initialize_parameterizer...\n");
    t1 = GetCurrentTime64();
    field.initialize_parameterizer(faces);
    t2 = GetCurrentTime64();
    printf("Use %lf seconds\n", (t2 - t1) * 1e-3);

    if (field.flag_preserve_boundary) {
        printf("Add boundary constrains...\n");
        Hierarchy &mRes = field.m_hierarchy;
        mRes.clearConstraints();
        for (uint32_t i = 0; i < 3 * mRes.mF.cols(); ++i) {
            if (mRes.mE2E[i] == -1) {
                uint32_t i0 = mRes.mF(i % 3, i / 3);
                uint32_t i1 = mRes.mF((i + 1) % 3, i / 3);
                Vector3d p0 = mRes.mV[0].col(i0), p1 = mRes.mV[0].col(i1);
                Vector3d edge = p1 - p0;
                if (edge.squaredNorm() > 0) {
                    edge.normalize();
                    mRes.mCO[0].col(i0) = p0;
                    mRes.mCO[0].col(i1) = p1;
                    mRes.mCQ[0].col(i0) = mRes.mCQ[0].col(i1) = edge;
                    mRes.mCQw[0][i0] = mRes.mCQw[0][i1] = mRes.mCOw[0][i0] = mRes.mCOw[0][i1] =
                            1.0;
                }
            }
        }
        mRes.propagateConstraints();
    }

    printf("Solve Orientation Field...\n");
    t1 = GetCurrentTime64();

    Optimizer::optimize_orientations(field.m_hierarchy);
    field.find_orientation_singularities();
    t2 = GetCurrentTime64();
    printf("Use %lf seconds\n", (t2 - t1) * 1e-3);

    if (field.flag_adaptive_scale == 1) {
        printf("Estimate Slop...\n");
        t1 = GetCurrentTime64();
        field.estimate_slope();
        t2 = GetCurrentTime64();
        printf("Use %lf seconds\n", (t2 - t1) * 1e-3);
    }
    printf("Solve for scale...\n");
    t1 = GetCurrentTime64();
    Optimizer::optimize_scale(field.m_hierarchy, field.rho, field.flag_adaptive_scale);
    field.flag_adaptive_scale = 1;
    t2 = GetCurrentTime64();
    printf("Use %lf seconds\n", (t2 - t1) * 1e-3);

    printf("Solve for position field...\n");
    t1 = GetCurrentTime64();
    Optimizer::optimize_positions(field.m_hierarchy, field.flag_adaptive_scale);

    field.find_position_singularities();
    t2 = GetCurrentTime64();
    printf("Use %lf seconds\n", (t2 - t1) * 1e-3);
    t1 = GetCurrentTime64();
    printf("Solve index map...\n");
    field.compute_index_map();
    t2 = GetCurrentTime64();
    printf("Indexmap Use %lf seconds\n", (t2 - t1) * 1e-3);
    printf("Writing the file...\n");


    if (output_obj.size() < 1) {
        assert(0);
        // field.save_to_obj((std::string(DATA_PATH) + "/result.obj").c_str());
    } else {
        field.save_to_obj(output_obj.c_str());
    }

//    const auto mesh_out = adapters::from_parametrizer_to_quad_mesh(field);
//    service.save_quadmesh_to_file(output_obj, mesh_out);
    printf("finish...\n");
    return 0;
}
