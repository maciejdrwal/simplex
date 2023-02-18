// Copyright (C) 2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 
#include "simplex.h"
#include "utils.h"
#include "logger.h"

#include <chrono>
#include <csignal>
#include <fstream>
#include <iostream>
#include <sstream>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#endif

#ifdef __linux__
#include <clbas.h>
#include <f77blas.h>
#endif

namespace simplex
{    
    int * LinearProgram::ipiv = NULL;

    std::chrono::time_point<std::chrono::steady_clock> t_start;

    // Print internal data structures on screen (for debugging).
    void LinearProgram::print_tableau() const
    {
        // Note that matrix_A is stored in a column-major order.
        std::stringstream ss;
        ss << "N=" << N << ", M=" << M << '\n';
        ss << "A=\n";
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                ss << matrix_A[j * M + i] << '\t';
            }
            ss << '\n';
        }
        ss << "b=\n";
        for (int i = 0; i < M; i++)
        {
            ss << vector_b[i] << '\t';
        }
        ss << "\nc=\n";
        for (int i = 0; i < N; i++)
        {
            ss << vector_c[i] << '\t';
        }
        ss << '\n';
        LOG(debug) << ss.str();
    }

    void print_matrix(double * matrix, int n_rows, int m_cols)
    {
        std::stringstream ss;
        for (int i = 0; i < n_rows; i++)
        {
            for (int j = 0; j < m_cols; j++)
            {
                ss << std::scientific << matrix[j * n_rows + i] << '\t';
            }
            ss << std::fixed << '\n';
            LOG(debug) << ss.str();
        }
    }

    template <typename T>
    void print_vector(T * vect, int n)
    {
        std::stringstream ss;
        for (int i = 0; i < n; i++)
        {
            ss << std::scientific << vect[i] << '\t';
        }
        ss << std::fixed << '\n';
        LOG(debug) << ss.str();
    }

    // Copies the data read by parser into internal data structures: A, b.
    // Variable names are mapped here to indices using two maps:
    // variable_name_to_id and variable_id_to_name.
    // This should be called just before simplex(), after presolve and all
    // other problem transformations.
    void LinearProgram::initialize_tableau()
    {
        variable_id_to_name.clear();
        variable_name_to_id.clear();
        int next_id = 0;
        int j = 0;

        for (const auto & constraint : constraints)
        {
            for (const auto & [name, coeff] : constraint.second.name_coeff)
            {
                if (variable_name_to_id.count(name) == 0)
                {
                    variable_id_to_name[next_id] = name;
                    variable_name_to_id[name] = next_id;
                    next_id++;
                }

                int i = variable_name_to_id[name] * M + j;
                matrix_A[i] = coeff;
            }
            vector_b[j++] = constraint.second.rhs;
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
                  << static_cast<long>(std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count()) << " ns)\n";
    }

    // The main procedure of the solver.
    // Initializes internal data structures and runs 2 phases of simplex algorithm.
    void LinearProgram::solve()
    {
        measure_time_start();

        M = constraints.size();
        N = var_lbnd.size();

        LOG(debug) << "problem size:\nN=" << N << ", M=" << M;
        LOG(debug) << "Allocating memory...";

        vector_b = std::make_unique<double[]>(M);
        vector_c = std::make_unique<double[]>(N + 2 * M);
        matrix_A = std::make_unique<double[]>(N * (N + 2 * M));

        std::fill_n(vector_b.get(), M, 0.0);
        std::fill_n(vector_c.get(), N + 2 * M, 0.0);
        std::fill_n(matrix_A.get(), N * (N + 2 * M), 0.0);

        std::set<std::string> init_basis; // Initial basis.

        // Add slack variables to inequality constraints, replacing them with equality.
        int num_of_slack_vars = 0;
        for (auto & constraint : constraints)
        {
            if (constraint.second.type != '=')
            {
                const std::string var_name("__SLACK" + utils::tostr<int>(num_of_slack_vars++));
                if (constraint.second.type == '<')
                {
                    constraint.second.name_coeff[var_name] = 1.0;
                }
                if (constraint.second.type == '>')
                {
                    constraint.second.name_coeff[var_name] = -1.0;
                }
                constraint.second.type = '=';
                add_variable(var_name);

                // If all constraints were inequalities then use slacks for initial basis.
                if (all_inequalities)
                {
                    init_basis.insert(var_name);
                }
                LOG(debug) << "added slack variable: " << var_name;
            }
        }
        int num_of_primal_vars = objective_name_coeff.size();

        N = var_lbnd.size(); // this might have been modified by adding slack variables

        if (!all_inequalities)
        {
            // Construct artificial variables for Phase I.
            int art_var_id = 0;
            LOG(debug) << "Not all constrains are inequalities A <= b, adding artificial variables.";
            for (auto & constraint : constraints)
            {
                if (constraint.second.rhs < 0.0)
                {
                    constraint.second.rhs *= -1.0;

                    for (auto & jt : constraint.second.name_coeff)
                    {
                        jt.second = -jt.second;
                    }
                    // Note: at this point there should be no inequality constraints
                    if (constraint.second.type == '<' || constraint.second.type == '>')
                    {
                        throw "All constraints must be equality at this point.";
                    }
                }
                const std::string var_name("__ARTIFICIAL" + utils::tostr<int>(art_var_id++));
                constraint.second.name_coeff[var_name] = 1.0;
                add_variable(var_name);
                init_basis.insert(var_name);
                LOG(debug) << "added artificial variable: " << var_name;
            }
            N = var_lbnd.size(); // this has changed after adding artificial variables

            // Solve Phase I LP.
            LOG(debug) << "*** PHASE I ***";
            initialize_tableau();
            // Objective function: sum of artificial variables.
            for (int i = 0; i < art_var_id; i++)
            {
                const std::string var_name("__ARTIFICIAL" + utils::tostr<int>(i));
                vector_c[variable_name_to_id[var_name]] = 1.0;
            }
            print_tableau();
            simplex(init_basis);

            // Remove artificial variables.
            int i = 0;
            for (auto & constraint : constraints)
            {
                const auto var_name = std::string("__ARTIFICIAL") + utils::tostr<int>(i++);
                int var_id = variable_name_to_id[var_name];

                // Check if some artificial variable remained in the basis.
                if (init_basis.count(var_name) > 0)
                {
                    // TODO: handle this case

                    // init_basis.erase(var_name);
                    // for (map<string, int>::iterator jt = variable_name_to_id.begin(); jt != variable_name_to_id.end(); ++jt) {
                    //     if (init_basis.count(jt->second) == 0 && (strncmp(jt->first.c_str(), "__ARTIFICIAL", 12) != 0)) {
                    //         init_basis.insert(jt->second);
                    //         printf("replaced variable %d (%s) by %d (%s) in the basis\n",  var_id, var_name.c_str(), jt->second, jt->first.c_str());
                    //         break;
                    //     }
                    // }

                    std::stringstream msg;
                    msg << "Artificial variable " << var_name << " in the basis.";
                    throw msg.str().c_str();
                }

                constraint.second.name_coeff.erase(var_name);
                variable_name_to_id.erase(var_name);
                variable_id_to_name.erase(var_id);
                remove_variable(var_name);
                LOG(debug) << "removed artificial variable: " << var_name;
            }
            N = var_lbnd.size(); // this has changed after removing artificial variables

            // TODO: ...
            std::fill_n(vector_b.get(), M, 0.0);
            std::fill_n(vector_c.get(), N + 2 * M, 0.0);
            std::fill_n(matrix_A.get(), N * (N + 2 * M), 0.0);
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

            for (const auto & [name, coeff] : objective_name_coeff)
            {
                int i = variable_name_to_id[name];
                vector_c[i] = (sense == 'm') ? coeff : -coeff;
            }

            // Solve Phase II LP.
            LOG(debug) << "*** PHASE II ***";
            print_tableau();
            simplex(init_basis);
        }

        print_elapsed_time();

        delete ipiv;
        ipiv = NULL;
    }

    // This functions solves a system of linear equations of the form:
    // A * x = b
    // The result is stored in _vector_bx, thus enough space must be allocated.
    // If factorize = true then also the following system is solved:
    // A^T * y = c
    // The result is stored in _vector_cy, thus enough space must be allocated.
    // If factorize = false then _matrix_A must contain LU factors of A on input.
    void LinearProgram::lineq_solve(double * _matrix_A, double * _vector_bx, double * _vector_cy, bool factorize = true)
    {
        int info;
        int one_int = 1;
        if (ipiv == nullptr) ipiv = new int[M];

        if (factorize)
        {
            // LU factorization: A = P * L * U
            // Arguments: num. of rows, num. of cols., matrix (on exit: LU factors), leading dimension of array, pivot indices, info.
            dgetrf_(&M, &M, _matrix_A, &M, ipiv, &info);
            if (info != 0)
            {
                // info < 0: argument "-info" has illegal value
                // info > 0: diagonal element at "info" of the factor U from the factorization is exactly 0
                LOG(debug) << "dgetrf_ failed with error code " << static_cast<int>(info);
                std::exit(1);
            }

            // Compute: A^T * y = c
            // Arguments: transpose, order of matrix, num. of r.h.sides, LU factors in single matrix (from dgetrf),
            // leading dimension of array, pivot indices from DGETRF, rhs (on exit: the solution), leading dimension of rhs, info.
            char transpose = 'T';
            dgetrs_(&transpose, &M, &one_int, _matrix_A, &M, ipiv, _vector_cy, &M, &info);
            if (info != 0)
            {
                // info < 0: argument "-info" has illegal value
                LOG(debug) << "dgetrs_ #1 failed with error code " << static_cast<int>(info);
                std::exit(1);
            }
        }

        // Compute: A * x = b
        char transpose = 'N';
        dgetrs_(&transpose, &M, &one_int, _matrix_A, &M, ipiv, _vector_bx, &M, &info);
        if (info != 0)
        {
            // info < 0: argument "-info" has illegal value
            LOG(debug) << "dgetrs_ #2 failed with error code " << static_cast<int>(info);
            std::exit(1);
        }
    }

    // Bland's rule for pivoting.
    // Choose entering index to be lexicographically first with s_j < 0.
    int LinearProgram::select_entering_variable_Bland(double * vector_c_N)
    {
        for (int i = 0; i < N - M; i++)
        {
            if (vector_c_N[i] < 0.0) return i;
        }
        return -1;
    }

    // Bland's rule for pivoting.
    // Choose leaving index to be lexicographically first with min{ x_i / d_i, d_i > 0 },
    // d = A_B^{-1} * A(entering_index)
    int LinearProgram::select_leaving_variable_Bland(double * vector_bx, double * vector_cy, double * matrix_A_B)
    {
        // Solve: A_B * d = A(entering_index). Note that matrix_A_B already contains LU factors.
        lineq_solve(matrix_A_B, vector_cy, nullptr, false);

        // cout << "solved d=\n";
        // print_vector<double>(vector_cy, M);

        int min_i = -1;
        double min_lbd = std::numeric_limits<double>::max();
        for (int i = 0; i < M; i++)
        {
            if (vector_cy[i] > 0.0)
            {
                double lbd = vector_bx[i] / vector_cy[i];
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
    int LinearProgram::select_entering_variable_most_neg(double * vector_c_N)
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
    int LinearProgram::select_leaving_variable_SUB(double * vector_bx, double * vector_cy, double * matrix_A_B, int entering_index)
    {
        // Solve: A_B * d = A(entering_index). Note that matrix_A_B already contains LU factors.
        // Result d is stored in vector_cy.
        lineq_solve(matrix_A_B, vector_cy, nullptr, false);

        bool need_substitution = false;
        int min_i = -1;
        double min_theta = std::numeric_limits<double>::max();
        for (int i = 0; i < M; i++)
        {
            if (vector_cy[i] > 0.0)
            {
                double t = vector_bx[i] / vector_cy[i];
                if (t < min_theta)
                {
                    min_i = i;
                    min_theta = t;
                }
            }
        }

        for (int i = 0; i < M; i++)
        {
            if (vector_cy[i] < 0.0)
            {
                double ub = var_ubnd[variable_id_to_name[basis[i]]];
                double t = (ub - vector_bx[i]) / (-vector_cy[i]);
                if (t < min_theta)
                {
                    min_i = i;
                    min_theta = t;
                    need_substitution = true;
                }
            }
        }

        double ub_s = var_ubnd[variable_id_to_name[non_basis[entering_index]]];

        if (ub_s < min_theta)
        {
            upper_bound_substitution(non_basis[entering_index], ub_s);
            return -2; // do not perform pivot
        }
        else
        {
            if (need_substitution)
            {
                int k = basis[min_i];
                upper_bound_substitution(k, var_ubnd[variable_id_to_name[k]]);
            }
        }

        return min_i; // perform usual pivot (except if min_i == -1)
    }

    void LinearProgram::upper_bound_substitution(int var_id, double ub)
    {
        // Modify the objective function: c_s x_s -> -c_s x_s
        vector_c[var_id] = -vector_c[var_id];

        // Modify constraints: a_{js} x_s -> -a_{js} x_s, b_j -> b_j - a_{js}*u_s
        for (int i = 0; i < M; i++)
        {
            double * a_ptr = &matrix_A[var_id * M + i];
            double a = *a_ptr;
            vector_b[i] = vector_b[i] - a * ub;
            *a_ptr = -a;
        }

        ub_substitutions[var_id] = !ub_substitutions[var_id];
    }

    void LinearProgram::solution_found(double * vector_bx, double * vector_cy, double * vector_c_B)
    {
        double cost = 0.0;
        solution.clear();

        for (int i = 0; i < basis.size(); i++)
        {
            int var_id = basis[i];
            const std::string & var_name = variable_id_to_name[var_id];

            // Reverse UB-substitution
            if (ub_substitutions[var_id])
            {
                LOG(debug) << "REVERSING UB-subst of " << var_name;
                upper_bound_substitution(var_id, var_ubnd[var_name]);
            }

            // Applying shifts
            for (auto kt = var_shifts.begin(); kt != var_shifts.end(); ++kt)
            {
                if (variable_name_to_id[kt->first] == var_id)
                {
                    LOG(debug) << "SHIFTING: basis var " << var_name << " by " << kt->second;
                    vector_bx[i] += kt->second;
                }
            }

            LOG(debug) << "var (basis) " << var_name << " = " << vector_bx[i];
            solution[var_name] = vector_bx[i];
            cost += vector_c_B[i] * vector_bx[i];
        }
        for (int i = 0; i < non_basis.size(); i++)
        {
            double q = 0.0;
            int var_id = non_basis[i];
            const std::string & var_name = variable_id_to_name[var_id];

            // Reverse UB-substitution
            if (ub_substitutions[var_id])
            {
                double ub = var_ubnd[var_name];
                upper_bound_substitution(var_id, ub);
                q += ub;
                LOG(debug) << "REVERSING UB-subst of " << var_name << " ub= " << ub;
            }

            // Applying shifts
            for (auto kt = var_shifts.begin(); kt != var_shifts.end(); ++kt)
            {
                if (variable_name_to_id[kt->first] == var_id)
                {
                    LOG(debug) << "SHIFTING: non-basis var " << var_name << " by " << kt->second;
                    q += kt->second;
                }
            }

            LOG(debug) << "var (non-basis) " << var_name << " = " << q;
            solution[var_name] = q;
            cost += vector_c[non_basis[i]] * q;
        }
        objective_value = (sense == 'm') ? cost : -cost;

        LOG(debug) << "OPTIMAL X FOUND.\nValue = " << objective_value;
    }

    // The (Revised) Simplex Algorithm.
    // arg_basis : on input: initial basis; on result: final basis
    int LinearProgram::simplex(std::set<std::string> & arg_basis)
    {
        double * matrix_A_B = new double[M * M];
        double * matrix_A_N = new double[M * (N - M)];
        double * vector_c_B = new double[M];
        double * vector_c_N = new double[N - M];
        double * vector_bx = new double[M];
        double * vector_cy = new double[N];

        basis.clear();
        non_basis.clear();
        for (int i = 0; i < N; i++)
        {
            if (arg_basis.count(variable_id_to_name[i]) == 0)
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
        unsigned long iteration_limit = 1000000000L;
        while (iteration_count < iteration_limit)
        {
            iteration_count++;
            LOG(debug) << "SIMPLEX iteration: " << iteration_count;
            LOG(debug) << "Basis: ";
            print_vector<int>(basis.data(), basis.size());

            // Split the matrix A into basis matrix (A_B) and non-basis matrix (A_N).
            // TODO: no need to copy whole vectors, as only 1 element has changed in basis.
            // However note that matrix_A, vector_b and vector_c might have changed
            // due to applying upper-bound substitutions.
            int jB = 0, jN = 0;
            for (auto it = basis.begin(); it != basis.end(); ++it)
            {
                std::copy_n(&matrix_A[(*it) * M], M, matrix_A_B + jB * M);
                vector_c_B[jB] = vector_c[*it];
                jB++;
            }
            for (auto it = non_basis.begin(); it != non_basis.end(); ++it)
            {
                std::copy_n(&matrix_A[(*it) * M], M, matrix_A_N + jN * M);
                vector_c_N[jN] = vector_c[*it];
                jN++;
            }

            LOG(debug) << "A_B=\n";
            print_matrix(matrix_A_B, M, M);
            LOG(debug) << "A_N=\n";
            print_matrix(matrix_A_N, M, N - M);
            LOG(debug) << "c_B=\n";
            print_vector<double>(vector_c_B, M);
            LOG(debug) << "c_N=\n";
            print_vector<double>(vector_c_N, N - M);

            // Compute x = A_B^{-1} * b and y = (A_B^T)^{-1} * c_B
            std::copy_n(vector_b.get(), M, vector_bx);
            std::copy_n(vector_c_B, M, vector_cy);
            lineq_solve(matrix_A_B, vector_bx, vector_cy);

            LOG(debug) << "solved x=\n";
            print_vector<double>(vector_bx, M);
            LOG(debug) << "solved y=\n";
            print_vector<double>(vector_cy, M);

            // Pricing (reduced costs): s = c_N - (A_N)^T * y
            // Multiply matrix by vector: y = alpha * A * x + beta * y.
            // Arguments: storage order, transpose, num. of rows, num. of cols., alpha, matrix A, l.d.a. of A, vect. x, incx, beta, y, incy.
            // NOTE 1: the first argument "storage order" is not in the original CBLAS specification
            // NOTE 2: the second argument should be 'T' according to the original CBLAS specification
            cblas_dgemv(CblasColMajor, CblasTrans, M, N - M, -1.0, matrix_A_N, M, vector_cy, 1, 1.0, vector_c_N, 1);

            // Now vector_c_N contains the result s.

            LOG(debug) << "reduced costs s=\n";
            print_vector<double>(vector_c_N, N - M);

            // 1) Select the entering variable.
            int entering_index = select_entering_variable_Bland(vector_c_N);
            // int entering_index = select_entering_variable_most_neg(vector_c_N);

            if (entering_index == -1)
            {
                solution_found(vector_bx, vector_cy, vector_c_B);
                arg_basis.clear();
                for (int i = 0; i < M; i++)
                {
                    arg_basis.insert(variable_id_to_name[basis[i]]);
                }
                break;
            }

            // Copy the entering column A_q into vector_cy.
            std::copy_n(&matrix_A[non_basis[entering_index] * M], M, vector_cy);

            // 2) Select the leaving variable.
            // int leaving_index = select_leaving_variable_Bland(vector_bx, vector_cy, matrix_A_B);
            int leaving_index = select_leaving_variable_SUB(vector_bx, vector_cy, matrix_A_B, entering_index);

            if (leaving_index == -1)
            {
                LOG(debug) << "LP IS UNBOUNDED.\n";
                break;
            }

            if (leaving_index != -2)
            {
                // Update basis.
                int tmp = basis[leaving_index];
                basis[leaving_index] = non_basis[entering_index];
                non_basis[entering_index] = tmp;
            }
        }

        if (iteration_count == iteration_limit)
        {
            throw "MAXIMUM NUMBER OF SIMPLEX ITERATIONS REACHED";
        }

        delete[] matrix_A_B;
        delete[] matrix_A_N;
        delete[] vector_c_B;
        delete[] vector_c_N;
        delete[] vector_bx;
        delete[] vector_cy;

        return 0;
    }

    // Write the solution to text file.
    void LinearProgram::write(const std::string & filename) const
    {
        LOG(debug) << "Writing solution...\n";

        std::ofstream outfile;
        outfile.open(filename);

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
    void LinearProgram::add_variable(const std::string & var_name)
    {
        var_lbnd[var_name] = 0.0; // set default bounds
        var_ubnd[var_name] = std::numeric_limits<double>::max();
    }

    void LinearProgram::remove_variable(const std::string & var_name)
    {
        var_lbnd.erase(var_name);
        var_ubnd.erase(var_name);
    }

    bool LinearProgram::has_variable(const std::string & var_name) const
    {
        return (0 != var_lbnd.count(var_name));
    }
}

