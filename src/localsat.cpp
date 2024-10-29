#include <deque>
#include <memory>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "localsat.h"
#include "dedge.h"
#include "field-math.h"

namespace qflow {

    const int max_depth = 0;

    using namespace Eigen;

    SolverStatus RunCNF(
            const std::string &fin_name,
            int n_variable,
            int timeout,
            const std::vector<std::vector<int>> &sat_clause,
            std::vector<int> &value
    ) {
        int n_sat_variable = 3 * n_variable;
        auto fout_name = fin_name + ".result.txt";

        FILE *fout = fopen(fin_name.c_str(), "w");
        fprintf(fout, "p cnf %d %d\n", n_sat_variable, (int) sat_clause.size());
        for (auto &c: sat_clause) {
            for (auto e: c) fprintf(fout, "%d ", e);
            fputs("0\n", fout);
        }
        fclose(fout);

        char cmd[100];
        snprintf(cmd, 99, "rm %s > /dev/null 2>&1", fout_name.c_str());
        system(cmd);
        snprintf(cmd, 99, "timeout %d minisat %s %s > /dev/null 2>&1", timeout, fin_name.c_str(),
                 fout_name.c_str());
        int exit_code = system(cmd);

        FILE *fin = fopen(fout_name.c_str(), "r");
        char buf[16] = {0};
        fscanf(fin, "%15s", buf);
//        lprintf("  MiniSAT:");
        if (strcmp(buf, "SAT") != 0) {
            fclose(fin);

            if (exit_code == 124) {
//                lprintf("       Timeout! ");
                return SolverStatus::Timeout;
            }
//            lprintf(" Unsatisfiable! ");
            return SolverStatus::Unsat;
        };

//        lprintf("   Satisfiable! ");
        for (int i = 0; i < n_variable; ++i) {
            int sign[3];
            fscanf(fin, "%d %d %d", sign + 0, sign + 1, sign + 2);

            int nvalue = -2;
            for (int j = 0; j < 3; ++j) {
                assert(abs(sign[j]) == 3 * i + j + 1);
                if ((sign[j] > 0) == (value[i] != j - 1)) {
                    assert(nvalue == -2);
                    nvalue = j - 1;
                }
            }
            value[i] = nvalue;
        }
        fclose(fin);

        return SolverStatus::Sat;
    }

    SolverStatus SolveSatProblem(
            int n_variable,
            std::vector<int> &value,
            const std::vector<bool> flexible,  // NOQA
            const std::vector<Vector3i> &variable_eq,
            const std::vector<Vector3i> &constant_eq,
            const std::vector<Vector4i> &variable_ge,
            const std::vector<Vector2i> &constant_ge,
            int timeout
    ) {
        for (int v: value) assert(-1 <= v && v <= +1);

        auto VAR = [&](int i, int v) {
            int index = 1 + 3 * i + v + 1;
            // We initialize the SAT problem by setting all the variable to false.
            // This is because minisat by default will try false first.
            if (v == value[i]) index = -index;
            return index;
        };

        int n_flexible = 0;
        std::vector<std::vector<int>> sat_clause;
        std::vector<bool> sat_ishard;

        auto add_clause = [&](const std::vector<int> &clause, bool hard) {
            sat_clause.push_back(clause);
            sat_ishard.push_back(hard);
        };

        for (int i = 0; i < n_variable; ++i) {
            add_clause({-VAR(i, -1), -VAR(i, 0)}, true);
            add_clause({-VAR(i, +1), -VAR(i, 0)}, true);
            add_clause({-VAR(i, -1), -VAR(i, +1)}, true);
            add_clause({VAR(i, -1), VAR(i, 0), VAR(i, +1)}, true);
            if (!flexible[i]) {
                add_clause({VAR(i, value[i])}, true);
            } else {
                ++n_flexible;
            }
        }

        for (int i = 0; i < (int) variable_eq.size(); ++i) {
            auto &var = variable_eq[i];
            auto &cst = constant_eq[i];
            for (int v0 = -1; v0 <= 1; ++v0)
                for (int v1 = -1; v1 <= 1; ++v1)
                    for (int v2 = -1; v2 <= 1; ++v2)
                        if (cst[0] * v0 + cst[1] * v1 + cst[2] * v2 != 0) {
                            add_clause({-VAR(var[0], v0), -VAR(var[1], v1), -VAR(var[2], v2)}, true);
                        }
        }

        for (int i = 0; i < (int) variable_ge.size(); ++i) {
            auto &var = variable_ge[i];
            auto &cst = constant_ge[i];
            for (int v0 = -1; v0 <= 1; ++v0)
                for (int v1 = -1; v1 <= 1; ++v1)
                    for (int v2 = -1; v2 <= 1; ++v2)
                        for (int v3 = -1; v3 <= 1; ++v3)
                            if (cst[0] * v0 * v1 - cst[1] * v2 * v3 < 0) {
                                add_clause({-VAR(var[0], v0), -VAR(var[1], v1), -VAR(var[2], v2),
                                            -VAR(var[3], v3)},
                                           false);
                            }
        }

        int nflip_before = 0, nflip_after = 0;
        for (int i = 0; i < (int) variable_ge.size(); ++i) {
            auto &var = variable_ge[i];
            auto &cst = constant_ge[i];
            if (value[var[0]] * value[var[1]] * cst[0] - value[var[2]] * value[var[3]] * cst[1] < 0)
                nflip_before++;
        }

//        lprintf("  [SAT] nvar: %6d nflip: %3d ", n_flexible * 2, nflip_before);
        auto rcnf = RunCNF("test.out", n_variable, timeout, sat_clause, value);

        for (int i = 0; i < (int) variable_eq.size(); ++i) {
            auto &var = variable_eq[i];
            auto &cst = constant_eq[i];
            assert(cst[0] * value[var[0]] + cst[1] * value[var[1]] + cst[2] * value[var[2]] == 0);
        }
        for (int i = 0; i < (int) variable_ge.size(); ++i) {
            auto &var = variable_ge[i];
            auto &cst = constant_ge[i];
            int area = value[var[0]] * value[var[1]] * cst[0] - value[var[2]] * value[var[3]] * cst[1];
            if (area < 0) ++nflip_after;
        }
//        lprintf("nflip: %3d\n", nflip_after);
        return rcnf;
    }
}
