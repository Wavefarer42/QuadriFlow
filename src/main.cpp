#include <cstdlib>
#include <chrono>

#include "field-math.h"
#include "optimizer.h"
#include "parametrizer.h"

using namespace qflow;

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
    if (input_obj.size() >= 1) {
        field.load_from_obj(input_obj.c_str());
    } else {
        assert(0);
        // field.load_from_obj((std::string(DATA_PATH) + "/fertility.obj").c_str());
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
    printf("finish...\n");
    //	field.LoopFace(2);
    return 0;
}
