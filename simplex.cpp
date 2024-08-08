// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#include "Simplex.h"
#include "Utils.h"
#include "Logger.h"

#include <chrono>
#include <csignal>
#include <fstream>
#include <sstream>

static const std::string ARTIFICIAL = "_ARTIFICIAL_";

namespace simplex
{
    std::chrono::time_point<std::chrono::steady_clock> t_start;

    // Print internal data structures on screen (for debugging).
    void LinearProgram::print_tableau() const
    {
        // Note that matrix_A is stored in a column-major order.
        std::stringstream ss;
        ss << "A=\n" << matrix_A << '\n';
        ss << "b=\n" << vector_b << '\n';
        ss << "c=\n" << vector_c << '\n';
        LOG(debug) << ss.str();
    }

    // Copies the data read by parser into internal data structures.
    // This should be called just before simplex(), after presolve and all
    // other problem transformations.
    void LinearProgram::initialize_tableau()
    {
        const auto M = constraints.size();
        const auto N = objective_coeff.size();

        LOG(debug) << "problem size: N=" << N << ", M=" << M;

        vector_b.resize(M);
        vector_c.resize(N);
        matrix_A.resize(M, N);

        vector_b.setZero();
        vector_c.setZero();
        matrix_A.setZero();

        int j = 0;
        for (const auto & [_, constraint] : constraints)
        {
            for (const auto & [name, a] : constraint.name_to_coeff)
            {
                const Eigen::Index i = variable_name_to_id[name];
                matrix_A(j, i) = a;
            }
            vector_b[j++] = constraint.rhs;
        }

        objective_value = 0.0;
    }

    inline void measure_time_start()
    {
        t_start = std::chrono::steady_clock::now();
    }

    inline void print_elapsed_time()
    {
        auto t_end = std::chrono::steady_clock::now();
        LOG(debug) << "Elapsed time: "
                   << static_cast<long>(std::chrono::duration_cast<std::chrono::seconds>(t_end - t_start).count()) << " s ("
                   << static_cast<long>(std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count()) << " ms, "
                   << static_cast<long>(std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count()) << " us, "
                   << static_cast<long>(std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count()) << " ns)";
    }

    bool LinearProgram::add_slack_variables_for_inequality_constraints(Basis & basis)
    {
        int num_of_slack_vars = 0;
        bool all_inequalities = true;
        for (auto & [_, constraint] : constraints)
        {
            if (constraint.type != '=')
            {
                const std::string var_name("__SLACK" + utils::to_str<int>(num_of_slack_vars++));
                if (constraint.type == '<')
                {
                    constraint.add_term(var_name, 1.0);
                }
                if (constraint.type == '>')
                {
                    constraint.add_term(var_name, -1.0);
                }
                constraint.type = '=';
                const auto var_id = add_variable(var_name);

                // If all constraints were inequalities then use slacks for initial basis.
                if (all_inequalities)
                {
                    basis.insert(var_id);
                }
                LOG(debug) << "added slack variable: " << var_name << " (" << var_id << ")";
            }
            else
            {
                all_inequalities = false;
                basis.clear();
            }
        }
        LOG(debug) << "added " << num_of_slack_vars << " slack variables, all_inequalities=" << all_inequalities;
        return all_inequalities;
    }

    void LinearProgram::add_artificial_variables_for_first_phase(Basis & basis)
    {
        int art_var_id = 0;
        LOG(debug) << "Not all constrains are inequalities A <= b, adding artificial variables.";
        for (auto & [_, constraint] : constraints)
        {
            // Note: at this point there should be no inequality constraints
            if (constraint.type == '<' || constraint.type == '>')
            {
                throw "All constraints must be equality at this point.";
            }

            if (constraint.rhs < 0.0)
            {
                constraint.rhs *= -1.0;

                for (auto & [_, a] : constraint.name_to_coeff)
                {
                    a *= -1;
                }
            }
            const std::string var_name(ARTIFICIAL + utils::to_str<int>(art_var_id++));
            constraint.add_term(var_name, 1.0);
            const auto var_id = add_variable(var_name);
            basis.insert(var_id);
            LOG(debug) << "added artificial variable: " << var_name;
        }
    }

    // The main procedure of the solver.
    // Initializes internal data structures and runs 2 phases of simplex algorithm.
    void LinearProgram::solve()
    {
        measure_time_start();

        Basis init_basis; // Initial basis.

        const bool all_inequalities = add_slack_variables_for_inequality_constraints(init_basis);

        // note: number of variables N might have been modified by adding slack variables

        if (!all_inequalities)
        {
            // Construct artificial variables for Phase I.
            add_artificial_variables_for_first_phase(init_basis);

            // note: number of variables N has changed after adding artificial variables

            // Solve Phase I LP.
            LOG(debug) << "*** PHASE I ***";
            initialize_tableau();
            // Objective function: sum of artificial variables.
            for (int i = 0; i < art_var_id; i++)
            {
                const std::string var_name(ARTIFICIAL + utils::to_str<int>(i));
                vector_c[variable_name_to_id[var_name]] = 1.0;
            }
            print_tableau();
            simplex(init_basis);

            // Remove artificial variables.
            int i = 0;
            for (auto & [_, constraint] : constraints)
            {
                const auto var_name = std::string(ARTIFICIAL) + utils::to_str<int>(i++);
                auto var_id = variable_name_to_id[var_name];

                // Check if some artificial variable remained in the basis.
                if (init_basis.count(var_id) > 0)
                {
                    // TODO: handle this case

                    // init_basis.erase(var_name);
                    // for (map<string, int>::iterator jt = variable_name_to_id.begin(); jt != variable_name_to_id.end(); ++jt) {
                    //     if (init_basis.count(jt->second) == 0 && (strncmp(jt->first.c_str(), ARTIFICIAL, 12) != 0)) {
                    //         init_basis.insert(jt->second);
                    //         printf("replaced variable %d (%s) by %d (%s) in the basis\n",  var_id, var_name.c_str(), jt->second, jt->first.c_str());
                    //         break;
                    //     }
                    // }

                    std::stringstream msg;
                    msg << "Artificial variable " << var_name << " in the basis.";
                    throw msg.str().c_str();
                }

                constraint.remove_term(var_name);
                remove_variable(var_id);
                variable_name_to_id.erase(var_name);
                LOG(debug) << "removed artificial variable: " << var_name;
                --N;
            }

            // note: number of variables N has changed after removing artificial variables

            vector_b.setZero();
            vector_c.setZero();
            matrix_A.setZero();
        }
        else
        {
            LOG(debug) << "Using slacks for initial basis, skipping Phase I.";
        }

        if (objective_value > 0.0)
        {
            LOG(debug) << "LP is infeasible. Aborting.";
        }
        else
        {
            initialize_tableau();

            // Prepare the original objective function.
            for (const auto &[i, coeff] : objective_coeff)
            {
                vector_c[i] = (sense == 'm') ? coeff : -coeff;
            }

            // Solve Phase II LP.
            LOG(debug) << "*** PHASE II ***";
            print_tableau();
            simplex(init_basis);
        }

        print_elapsed_time();
    }

    // This functions solves a system of linear equations of the form:
    // A * x = b
    // The result is stored in _vector_bx, thus enough space must be allocated.
    // If factorize = true then also the following system is solved:
    // A^T * y = c
    // The result is stored in _vector_cy, thus enough space must be allocated.
    // If factorize = false then _matrix_A must contain LU factors of A on input.
    void LinearProgram::lineq_solve(Eigen::MatrixXd & matrix_A, Eigen::VectorXd & vector_bx, Eigen::VectorXd & vector_cy, bool factorize = true)
    {
        // int info;
        // int one_int = 1;
        // if (ipiv == nullptr) ipiv = new int[M];

        // if (factorize)
        //{
        //     // LU factorization: A = P * L * U
        //     // Arguments: num. of rows, num. of cols., matrix (on exit: LU factors), leading dimension of array, pivot indices, info.
        //     dgetrf_(&M, &M, _matrix_A, &M, ipiv, &info);
        //     if (info != 0)
        //     {
        //         // info < 0: argument "-info" has illegal value
        //         // info > 0: diagonal element at "info" of the factor U from the factorization is exactly 0
        //         LOG(debug) << "dgetrf_ failed with error code " << static_cast<int>(info);
        //         std::exit(1);
        //     }

        //    // Compute: A^T * y = c
        //    // Arguments: transpose, order of matrix, num. of r.h.sides, LU factors in single matrix (from dgetrf),
        //    // leading dimension of array, pivot indices from DGETRF, rhs (on exit: the solution), leading dimension of rhs, info.
        //    char transpose = 'T';
        //    dgetrs_(&transpose, &M, &one_int, _matrix_A, &M, ipiv, _vector_cy, &M, &info);
        //    if (info != 0)
        //    {
        //        // info < 0: argument "-info" has illegal value
        //        LOG(debug) << "dgetrs_ #1 failed with error code " << static_cast<int>(info);
        //        std::exit(1);
        //    }
        //}

        //// Compute: A * x = b
        // char transpose = 'N';
        // dgetrs_(&transpose, &M, &one_int, _matrix_A, &M, ipiv, _vector_bx, &M, &info);
        // if (info != 0)
        //{
        //     // info < 0: argument "-info" has illegal value
        //     LOG(debug) << "dgetrs_ #2 failed with error code " << static_cast<int>(info);
        //     std::exit(1);
        // }

        if (factorize)
        {
            const Eigen::VectorXd y = matrix_A.transpose().colPivHouseholderQr().solve(vector_cy);
            vector_cy = y;
        }

        const Eigen::VectorXd x = matrix_A.colPivHouseholderQr().solve(vector_bx);
        vector_bx = x;
    }

    // Bland's rule for pivoting.
    // Choose entering index to be lexicographically first with s_j < 0.
    int LinearProgram::select_entering_variable_Bland(const Eigen::VectorXd &vector_c_N) const
    {
        for (int i = 0; i < N - M; i++)
        {
            if (vector_c_N[i] < 0.0)
            {
                return i;
            }
        }
        return -1;
    }

    // Bland's rule for pivoting.
    // Choose leaving index to be lexicographically first with min{ x_i / d_i, d_i > 0 },
    // d = A_B^{-1} * A(entering_index)
    int LinearProgram::select_leaving_variable_Bland(Eigen::VectorXd &vector_bx, Eigen::VectorXd &vector_cy, Eigen::MatrixXd &matrix_A_B)
    {
        // Solve: A_B * d = A(entering_index). Note that matrix_A_B already contains LU factors.
        Eigen::VectorXd unused;
        lineq_solve(matrix_A_B, vector_cy, unused, false);

        int min_i = -1;
        double min_lbd = std::numeric_limits<double>::max();
        for (int i = 0; i < M; i++)
        {
            if (vector_cy[i] > 0.0)
            {
                const double lbd = vector_bx[i] / vector_cy[i];
                if (lbd < min_lbd)
                {
                    min_i = i;
                    min_lbd = lbd;
                }
            }
        }
        return min_i;
    }

    // Choose the most negative value among s_j < 0.
    int LinearProgram::select_entering_variable_most_neg(const Eigen::VectorXd &vector_c_N) const
    {
        int min_i = -1;
        double min_value = std::numeric_limits<double>::max();
        for (int i = 0; i < N - M; i++)
        {
            if (vector_c_N[i] < 0.0 && vector_c_N[i] < min_value)
            {
                min_i = i;
                min_value = vector_c_N[i];
            }
        }
        return min_i;
    }

    // Simple Upper Bound rule for pivoting.
    int LinearProgram::select_leaving_variable_SUB(double *vector_bx, double *vector_cy, double *matrix_A_B, int entering_index)
    {
        // Solve: A_B * d = A(entering_index). Note that matrix_A_B already contains LU factors.
        // Result d is stored in vector_cy.
        // lineq_solve(matrix_A_B, vector_cy, nullptr, false);

        bool need_substitution = false;
        int min_i = -1;
        // double min_theta = std::numeric_limits<double>::max();
        // for (int i = 0; i < M; i++)
        //{
        //     if (vector_cy[i] > 0.0)
        //     {
        //         double t = vector_bx[i] / vector_cy[i];
        //         if (t < min_theta)
        //         {
        //             min_i = i;
        //             min_theta = t;
        //         }
        //     }
        // }

        // for (int i = 0; i < M; i++)
        //{
        //     if (vector_cy[i] < 0.0)
        //     {
        //         double ub = var_ubnd[variable_id_to_name[basis[i]]];
        //         double t = (ub - vector_bx[i]) / (-vector_cy[i]);
        //         if (t < min_theta)
        //         {
        //             min_i = i;
        //             min_theta = t;
        //             need_substitution = true;
        //         }
        //     }
        // }

        // double ub_s = var_ubnd[variable_id_to_name[non_basis[entering_index]]];

        // if (ub_s < min_theta)
        //{
        //     upper_bound_substitution(non_basis[entering_index], ub_s);
        //     return -2; // do not perform pivot
        // }
        // else
        //{
        //     if (need_substitution)
        //     {
        //         int k = basis[min_i];
        //         upper_bound_substitution(k, var_ubnd[variable_id_to_name[k]]);
        //     }
        // }

        return min_i; // perform usual pivot (except if min_i == -1)
    }

    void LinearProgram::upper_bound_substitution(Eigen::Index var_id, double ub)
    {
        // Modify the objective function: c_s x_s -> -c_s x_s
        vector_c[var_id] = -vector_c[var_id];

        // Modify constraints: a_{js} x_s -> -a_{js} x_s, b_j -> b_j - a_{js}*u_s
        for (auto i = 0u; i < M; i++)
        {
            const auto a = matrix_A(i, var_id);
            vector_b[i] = vector_b[i] - a * ub;
            matrix_A(i, var_id) = -a;
        }

        ub_substitutions[var_id] = !ub_substitutions[var_id];
    }

    void LinearProgram::solution_found(Eigen::VectorXd &vector_bx, const Eigen::VectorXd &vector_c_B)
    {
        double cost = 0.0;
        solution.clear();

        for (auto i = 0u; i < basis.size(); ++i)
        {
            const auto var_id = basis[i];
            const std::string &var_name = variable_id_to_name[var_id];

            // Reverse UB-substitution
            if (ub_substitutions[var_id])
            {
                LOG(debug) << "REVERSING UB-subst of " << var_name;
                upper_bound_substitution(var_id, var_ubnd[var_id]);
            }

            // Applying shifts
            for (auto kt = var_shifts.begin(); kt != var_shifts.end(); ++kt)
            {
                if (kt->first == var_id)
                {
                    LOG(debug) << "SHIFTING: basis var " << var_name << " by " << kt->second;
                    vector_bx[var_id] += kt->second;
                }
            }

            LOG(debug) << "var (basis) " << var_name << " = " << vector_bx[i];
            solution[var_name] = vector_bx[i];
            cost += vector_c_B[i] * vector_bx[i];
        }
        for (auto i = 0u; i < non_basis.size(); ++i)
        {
            const auto var_id = non_basis[i];
            double q = 0.0;
            const std::string &var_name = variable_id_to_name[var_id];

            // Reverse UB-substitution
            if (ub_substitutions[var_id])
            {
                double ub = var_ubnd[var_id];
                upper_bound_substitution(var_id, ub);
                q += ub;
                LOG(debug) << "REVERSING UB-subst of " << var_name << " ub= " << ub;
            }

            // Applying shifts
            for (auto kt = var_shifts.begin(); kt != var_shifts.end(); ++kt)
            {
                if (kt->first == var_id)
                {
                    LOG(debug) << "SHIFTING: non-basis var " << var_name << " by " << kt->second;
                    q += kt->second;
                }
            }

            LOG(debug) << "var (non-basis) " << var_name << " = " << q;
            solution[var_name] = q;
            cost += vector_c[var_id] * q;
        }
        objective_value = (sense == 'm') ? cost : -cost;

        LOG(debug) << "OPTIMAL X FOUND.\nValue = " << objective_value;
    }

    std::string print_basis(const std::vector<Eigen::Index> &basis)
    {
        std::stringstream ss;
        for (auto i : basis)
        {
            ss << i << ' ';
        }
        return ss.str();
    }

    // The (Revised) Simplex Algorithm.
    // arg_basis : on input: initial basis; on result: final basis
    int LinearProgram::simplex(std::set<Eigen::Index> &arg_basis)
    {
        std::stringstream ss;
        for (const auto var : arg_basis)
        {
            ss << var << ' ';
        }
        LOG(info) << "simplex basis: " << ss.str();

        Eigen::MatrixXd matrix_A_B(M, M);
        Eigen::MatrixXd matrix_A_N(M, N - M);
        Eigen::VectorXd vector_c_B(M);
        Eigen::VectorXd vector_c_N(N - M);

        basis.clear();
        non_basis.clear();
        for (auto i = 0; i < N; i++)
        {
            if (arg_basis.count(i) == 0)
            {
                non_basis.push_back(i);
            }
            else
            {
                basis.push_back(i);
            }
        }

        // Initialize upper-bound substitutions.
        ub_substitutions.assign(N, false);

        unsigned long iteration_count = 0L;
        constexpr unsigned long iteration_limit = 10L; // 1000000000L;
        while (iteration_count < iteration_limit)
        {
            iteration_count++;
            LOG(debug) << "SIMPLEX iteration: " << iteration_count;
            LOG(debug) << "Basis: " << print_basis(basis);

            // Split the matrix A into basis matrix (A_B) and non-basis matrix (A_N).
            // TODO: no need to copy whole vectors, as only 1 element has changed in basis.
            // However note that matrix_A, vector_b and vector_c might have changed
            // due to applying upper-bound substitutions.
            int jB = 0, jN = 0;
            for (auto index : basis)
            {
                matrix_A_B.col(jB) = matrix_A.col(index);
                vector_c_B[jB] = vector_c[index];
                jB++;
            }
            for (auto index : non_basis)
            {
                matrix_A_N.col(jN) = matrix_A.col(index);
                vector_c_N[jN] = vector_c[index];
                jN++;
            }

            LOG(debug) << "A_B=\n" << matrix_A_B;
            LOG(debug) << "A_N=\n" << matrix_A_N;
            LOG(debug) << "c_B=\n" << vector_c_B;
            LOG(debug) << "c_N=\n" << vector_c_N;

            // Compute x = A_B^{-1} * b and y = (A_B^T)^{-1} * c_B
            Eigen::VectorXd vector_bx = vector_b;
            Eigen::VectorXd vector_cy = vector_c_B;
            lineq_solve(matrix_A_B, vector_bx, vector_cy);

            LOG(debug) << "solved x=\n" << vector_bx;
            LOG(debug) << "solved y=\n" << vector_cy;

            // Pricing (reduced costs): s = c_N - (A_N)^T * y
            // Multiply matrix by vector: y = alpha * A * x + beta * y.
            // Arguments: storage order, transpose, num. of rows, num. of cols., alpha, matrix A, l.d.a. of A, vect. x, incx, beta, y, incy.
            // NOTE 1: the first argument "storage order" is not in the original CBLAS specification
            // NOTE 2: the second argument should be 'T' according to the original CBLAS
            // cblas_dgemv(CblasColMajor, CblasTrans, M, N - M, -1.0, matrix_A_N.get(), M, vector_cy.get(), 1, 1.0, vector_c_N.get(), 1);
            // Now vector_c_N contains the result s.

            Eigen::VectorXd s = vector_c_N - matrix_A_N.transpose() * vector_cy;

            LOG(debug) << "reduced costs s=\n"
                       << s;

            // 1) Select the entering variable.
            const int entering_index = select_entering_variable_Bland(s);
            // int entering_index = select_entering_variable_most_neg(s);

            if (entering_index == -1)
            {
                solution_found(vector_bx, vector_c_B);
                arg_basis.clear();
                for (auto i = 0; i < M; i++)
                {
                    arg_basis.insert(basis[i]);
                }
                break;
            }

            // Copy the entering column A_q into vector_cy.
            // std::copy_n(&matrix_A[non_basis[entering_index] * M], M, vector_cy.get());
            vector_cy = matrix_A.col(non_basis[entering_index]);

            // 2) Select the leaving variable.
            const auto leaving_index = select_leaving_variable_Bland(vector_bx, vector_cy, matrix_A_B);
            // int leaving_index = select_leaving_variable_SUB(vector_bx, vector_cy, matrix_A_B, entering_index);

            if (leaving_index == -1)
            {
                LOG(debug) << "LP IS UNBOUNDED.\n";
                break;
            }

            if (leaving_index != -2)
            {
                // Update basis.
                const auto tmp = basis[leaving_index];
                basis[leaving_index] = non_basis[entering_index];
                non_basis[entering_index] = tmp;
            }
        }

        if (iteration_count == iteration_limit)
        {
            throw "MAXIMUM NUMBER OF SIMPLEX ITERATIONS REACHED";
        }

        return 0;
    }

    // Write the solution to text file.
    void LinearProgram::write(const std::string &filename) const
    {
        LOG(debug) << "Writing solution...\n";

        std::ofstream outfile;
        outfile.open(filename);

        outfile << std::fixed;
        outfile << "<?xml version = \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << '\n';
        outfile << "<TestSolution>\n<header objectiveValue=\"" << objective_value << "\"/>" << '\n';

        for (auto it = solution.begin(); it != solution.end(); it++)
        {
            outfile << "<variable name=\"" << it->first.c_str() << "\" value=\"" << it->second << "\"/>" << '\n';
        }
        outfile << "</TestSolution>" << '\n';
        outfile.close();
    }

    // Note: each variable has ID that corresponds to that variable's position
    // in the computation array vector_c, and is also used, e.g., to
    // indicate basic/non-basic variables. Thus it is currently not allowed
    // to add variables after removing some other variables.
    // TODO: change this in order to allow problem manipulation after solving.
    ///@param var_name - name of the new variable
    ///@param coeff - the coefficient of the variable in the objective function
    ///@return the index of the new variable
    Eigen::Index LinearProgram::add_variable(const std::string & var_name, double coeff)
    {
        const auto var_id = next_variable_id++;
        variable_name_to_id[var_name] = var_id;
        variable_id_to_name[var_id] = var_name;
        objective_coeff[var_id] = coeff;
        var_lbnd[var_id] = 0.0; // set default bounds
        var_ubnd[var_id] = std::numeric_limits<double>::max();

        return var_id;
    }

    void LinearProgram::remove_variable(Eigen::Index var_id)
    {
        objective_coeff.erase(var_id);
        var_lbnd.erase(var_id);
        var_ubnd.erase(var_id);
        variable_id_to_name.erase(var_id);
    }

    bool LinearProgram::has_variable(const std::string &var_name) const
    {
        return variable_name_to_id.count(var_name) > 0;
    }
}
