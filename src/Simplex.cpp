// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#include "Simplex.h"
#include "Logger.h"

#include <chrono>
#include <csignal>
#include <fstream>
#include <sstream>

namespace
{
    constexpr double ALMOST_ZERO = 1e-10;

    std::chrono::time_point<std::chrono::steady_clock> t_start;

    void measure_time_start() { t_start = std::chrono::steady_clock::now(); }

    void print_elapsed_time()
    {
        auto t_end = std::chrono::steady_clock::now();
        LOG(debug) << "Elapsed time: "
                   << static_cast<long>(std::chrono::duration_cast<std::chrono::seconds>(t_end - t_start).count())
                   << " s ("
                   << static_cast<long>(std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count())
                   << " ms, "
                   << static_cast<long>(std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count())
                   << " us, "
                   << static_cast<long>(std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count())
                   << " ns)";
    }
}

namespace simplex
{
    // The main procedure of the solver.
    // Initializes internal data structures and runs 2 phases of simplex algorithm.
    void Simplex::solve()
    {
        measure_time_start();

        Basis init_basis;  // Initial basis.

        const bool all_inequalities = m_lp.add_slack_variables_for_inequality_constraints(init_basis);

        // note: number of variables N might have been modified by adding slack variables

        if (!all_inequalities)
        {
            // Construct artificial variables for Phase I.
            const auto art_var_id = m_lp.add_artificial_variables_for_first_phase(init_basis);

            // note: number of variables N has changed after adding artificial variables

            // Solve Phase I LP.
            LOG(debug) << "*** PHASE I ***";
            m_lp.initialize_tableau();

            // Objective function: sum of artificial variables.
            m_lp.vector_c.setZero();
            for (int i = 0; i < art_var_id; i++)
            {
                const auto var_name = LinearProgram::get_artificial_variable(i);
                m_lp.vector_c[m_lp.variable_name_to_id[var_name]] = 1.0;
            }
            m_lp.print_tableau();
            simplex(init_basis);

            // Remove artificial variables.
            int i = 0;
            for (auto & [_, constraint] : m_lp.constraints)
            {
                const auto var_name = LinearProgram::get_artificial_variable(i++);
                auto var_id = m_lp.variable_name_to_id[var_name];

                // Check if some artificial variable remained in the basis.
                if (std::find(init_basis.begin(), init_basis.end(), var_id) != init_basis.end())
                {
                    // TODO: handle this case

                    // init_basis.erase(var_name);
                    // for (map<string, int>::iterator jt = variable_name_to_id.begin(); jt !=
                    // variable_name_to_id.end(); ++jt) {
                    //     if (init_basis.count(jt->second) == 0 && (strncmp(jt->first.c_str(), ARTIFICIAL, 12) != 0)) {
                    //         init_basis.insert(jt->second);
                    //         printf("replaced variable %d (%s) by %d (%s) in the basis\n",  var_id, var_name.c_str(),
                    //         jt->second, jt->first.c_str()); break;
                    //     }
                    // }

                    std::stringstream msg;
                    msg << "Artificial variable " << var_name << " in the basis.";
                    throw msg.str().c_str();
                }

                constraint.remove_term(var_name);
                m_lp.remove_variable(var_id);
                m_lp.variable_name_to_id.erase(var_name);
                LOG(debug) << "removed artificial variable: " << var_name;
            }

            // note: number of variables N has changed after removing artificial variables

            m_lp.vector_b.setZero();
            m_lp.vector_c.setZero();
            m_lp.matrix_A.setZero();
        }
        else
        {
            LOG(debug) << "Using slacks for initial basis, skipping Phase I.";
        }

        if (m_objective_value > 0.0)
        {
            LOG(debug) << "LP is infeasible. Aborting.";
        }
        else
        {
            m_lp.initialize_tableau();

            // Solve Phase II LP.
            LOG(debug) << "*** PHASE II ***";
            m_lp.print_tableau();
            simplex(init_basis);
        }

        print_elapsed_time();
    }

    // Bland's rule for pivoting.
    // Choose entering index to be lexicographically first with s_j < 0.
    int Simplex::select_entering_variable_Bland(const Eigen::VectorXd & s) const
    {
        const auto M = m_lp.constraints.size();
        const auto N = m_lp.objective_coeff.size();
        for (auto i = 0u; i < N - M; i++)
        {
            if (s[i] < -ALMOST_ZERO)
            {
                return i;
            }
        }
        return -1;
    }

    // Bland's rule for pivoting.
    // Choose leaving index to be lexicographically first with min{ x_i / d_i, d_i > 0 },
    // d = A_B^{-1} * A(entering_index)
    int Simplex::select_leaving_variable_Bland(const Eigen::VectorXd & x, const Eigen::VectorXd & d) const
    {
        const auto M = m_lp.constraints.size();
        int min_i = -1;
        double min_lbd = std::numeric_limits<double>::max();
        for (auto i = 0u; i < M; i++)
        {
            if (d[i] > ALMOST_ZERO)
            {
                const double lbd = x[i] / d[i];
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
    int Simplex::select_entering_variable_most_neg(const Eigen::VectorXd & s) const
    {
        int min_i = -1;
        double min_value = std::numeric_limits<double>::max();
        const auto M = m_lp.constraints.size();
        const auto N = m_lp.objective_coeff.size();
        for (auto i = 0u; i < N - M; i++)
        {
            if (s[i] < -ALMOST_ZERO && s[i] < min_value)
            {
                min_i = i;
                min_value = s[i];
            }
        }
        return min_i;
    }

    int Simplex::select_leaving_variable(const Eigen::VectorXd & x, const Eigen::VectorXd & d,
                                         const Basis & basis, const Basis & non_basis,
                                         int entering_index) const
    {
        if (m_lp.has_non_trivial_upper_bounds())
        {
            return select_leaving_variable_SUB(x, d, basis, non_basis, entering_index);
        }
        return select_leaving_variable_Bland(x, d);
    }

    // Simple Upper Bound rule for pivoting.
    int Simplex::select_leaving_variable_SUB(const Eigen::VectorXd & x, const Eigen::VectorXd & d, const Basis & basis, const Basis & non_basis,
                                             int entering_index) const
    {
        bool need_substitution = false;
        int min_i = -1;
        double min_theta = std::numeric_limits<double>::max();
        const auto M = m_lp.constraints.size();
        for (auto i = 0u; i < M; i++)
        {
            if (d[i] > ALMOST_ZERO)
            {
                double t = x[i] / d[i];
                if (t < min_theta)
                {
                    min_i = i;
                    min_theta = t;
                }
            }
        }

        for (auto i = 0u; i < M; i++)
        {
            if (d[i] < -ALMOST_ZERO)
            {
                const auto it = m_lp.var_ubnd.find(basis[i]);
                if (it == m_lp.var_ubnd.end())
                {
                    throw "No upper bound found for basis element: " + std::to_string(basis[i]);
                }

                double ub = it->second;
                double t = (ub - x[i]) / (-d[i]);
                if (t < min_theta)
                {
                    min_i = i;
                    min_theta = t;
                    need_substitution = true;
                }
            }
        }

        const auto it = m_lp.var_ubnd.find(non_basis[entering_index]);
        double ub_s = it->second;

        if (ub_s < min_theta)
        {
            m_lp.upper_bound_substitution(non_basis[entering_index], ub_s);
            return -2;  // do not perform pivot
        }
        if (need_substitution)
        {
            const auto k = basis[min_i];

            const auto it1 = m_lp.var_ubnd.find(k);
            if (it1 == m_lp.var_ubnd.end())
            {
                throw "No upper bound found for basis element: " + std::to_string(k);
            }

            m_lp.upper_bound_substitution(k, it1->second);
        }

        return min_i;  // perform usual pivot (except if min_i == -1)
    }

    void Simplex::solution_found(Eigen::VectorXd vector_bx, const Eigen::VectorXd & vector_c_B, const Basis & basis, const Basis & non_basis)
    {
        double cost = 0.0;
        m_solution.clear();

        for (size_t i = 0u; i < basis.size(); ++i)
        {
            const auto var_id = basis[i];
            const std::string & var_name = m_lp.variable_id_to_name[var_id];

            // Reverse UB-substitution
            if (m_lp.ub_substitutions[var_id])
            {
                LOG(debug) << "REVERSING UB-subst of " << var_name;
                m_lp.upper_bound_substitution(var_id, m_lp.var_ubnd[var_id]);
            }

            // Applying shifts

            auto kt = m_lp.var_shifts.find(var_id);
            if (kt != m_lp.var_shifts.end())
            {
                LOG(debug) << "SHIFTING: basis var " << var_name << " by " << kt->second;
                vector_bx[i] += kt->second;
            }

            LOG(debug) << "var (basis) " << var_name << " = " << vector_bx[i];
            m_solution[var_name] = vector_bx[i];
            cost += vector_c_B[i] * vector_bx[i];
        }

        for (auto var_id : non_basis)
        {
            double q = 0.0;
            const std::string & var_name = m_lp.variable_id_to_name[var_id];

            // Reverse UB-substitution
            if (m_lp.ub_substitutions[var_id])
            {
                double ub = m_lp.var_ubnd[var_id];
                m_lp.upper_bound_substitution(var_id, ub);
                q += ub;
                LOG(debug) << "REVERSING UB-subst of " << var_name << " ub= " << ub;
            }

            // Applying shifts
            for (auto kt = m_lp.var_shifts.begin(); kt != m_lp.var_shifts.end(); ++kt)
            {
                if (kt->first == var_id)
                {
                    LOG(debug) << "SHIFTING: non-basis var " << var_name << " by " << kt->second;
                    q += kt->second;
                }
            }

            LOG(debug) << "var (non-basis) " << var_name << " = " << q;
            m_solution[var_name] = q;
            cost += m_lp.vector_c[var_id] * q;
        }
        m_objective_value = (m_lp.sense == 'm') ? cost : -cost;

        LOG(debug) << "OPTIMAL X FOUND.\nValue = " << m_objective_value;
    }

    std::string print_basis(const std::vector<Eigen::Index> & basis)
    {
        std::stringstream ss;
        for (auto i : basis)
        {
            ss << i << ' ';
        }
        return ss.str();
    }

    int Simplex::select_entering_variable(const Eigen::VectorXd & s) const
    {
        return select_entering_variable_Bland(s);
        // return select_entering_variable_most_neg(s);
    }

    // The (Revised) Simplex Algorithm.
    // arg_basis : on input: initial basis; on result: final basis
    int Simplex::simplex(Basis & arg_basis)
    {
        const auto M = m_lp.constraints.size();
        const auto N = m_lp.objective_coeff.size();

        Eigen::MatrixXd matrix_A_B(M, M);
        Eigen::MatrixXd matrix_A_N(M, N - M);
        Eigen::VectorXd vector_c_B(M);
        Eigen::VectorXd vector_c_N(N - M);

        Basis basis;
        Basis non_basis;
        for (auto i = 0u; i < N; i++)
        {
            if (std::find(arg_basis.begin(), arg_basis.end(), i) != arg_basis.end())
            {
                basis.push_back(i);
            }
            else
            {
                non_basis.push_back(i);
            }
        }

        // Initialize upper-bound substitutions.
        m_lp.ub_substitutions.assign(N, false);

        unsigned long iteration_count = 0L;
        constexpr unsigned long iteration_limit = 100L;  // 1000000000L;
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
                matrix_A_B.col(jB) = m_lp.matrix_A.col(index);
                vector_c_B[jB] = m_lp.vector_c[index];
                jB++;
            }
            for (auto index : non_basis)
            {
                matrix_A_N.col(jN) = m_lp.matrix_A.col(index);
                vector_c_N[jN] = m_lp.vector_c[index];
                jN++;
            }

            LOG(debug) << "A_B=\n" << matrix_A_B;
            LOG(debug) << "A_N=\n" << matrix_A_N;
            LOG(debug) << "c_B=\n" << vector_c_B;
            LOG(debug) << "c_N=\n" << vector_c_N;

            // Compute x = A_B^{-1} * b and y = (A_B^T)^{-1} * c_B
            const Eigen::VectorXd x = matrix_A_B.colPivHouseholderQr().solve(m_lp.vector_b);
            const Eigen::VectorXd y = matrix_A_B.transpose().colPivHouseholderQr().solve(vector_c_B);

            LOG(debug) << "solved x=\n" << x;
            LOG(debug) << "solved y=\n" << y;

            // Pricing (reduced costs): s = c_N - (A_N)^T * y
            Eigen::VectorXd s = vector_c_N - matrix_A_N.transpose() * y;

            LOG(debug) << "reduced costs s=\n" << s;

            // 1) Select the entering variable.
            const int entering_index = select_entering_variable(s);
            if (entering_index == -1)
            {
                solution_found(x, vector_c_B, basis, non_basis);
                arg_basis = basis;
                break;
            }

            // Copy the entering column A_q into vector_cy.
            Eigen::VectorXd A_entering = m_lp.matrix_A.col(non_basis[entering_index]);

            // 2) Select the leaving variable.

            // Solve: A_B * d = A(entering_index). Note that matrix_A_B already contains LU factors.
            const Eigen::VectorXd d = matrix_A_B.colPivHouseholderQr().solve(A_entering);

            const auto leaving_index = select_leaving_variable(x, d, basis, non_basis, entering_index);
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

    // Write the m_solution to text file.
    void Simplex::write(const std::string & filename) const
    {
        LOG(debug) << "Writing solution...\n";

        std::ofstream outfile;
        outfile.open(filename);

        outfile << std::fixed;
        outfile << "<?xml version = \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << '\n';
        outfile << "<TestSolution>\n<header objectiveValue=\"" << m_objective_value << "\"/>" << '\n';

        for (const auto & [var_name, value] : m_solution)
        {
            outfile << "<variable name=\"" << var_name << "\" value=\"" << value << "\"/>" << '\n';
        }
        outfile << "</TestSolution>" << '\n';
        outfile.close();
    }
}  // namespace simplex
