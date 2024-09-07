// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#include "Simplex.h"
#include "Logger.h"
#include "Utils.h"

#include <csignal>
#include <fstream>
#include <sstream>

#include "InitialBasis.h"

namespace
{
    constexpr double ALMOST_ZERO = 1e-7;
}

namespace simplex
{
    // The main procedure of the solver.
    // Initializes internal data structures and runs 2 phases of simplex algorithm.
    void Simplex::solve()
    {
        MEASURE_TIME_START(Solver);

        InitialBasis basis_initializer(m_lp, *this);

        Basis init_basis = basis_initializer.get_basis();

        if (m_objective_value > 0.0)
        {
            LOG(debug) << "LP is infeasible. Aborting.";
        }
        else
        {
            m_lp.initialize_tableau();

            // Solve Phase II LP.
            LOG(debug) << "*** PHASE II ***";

            //m_lp.print_tableau();
            simplex(init_basis);
        }

        LOG(info) << PRINT_ELAPSED_TIME(Solver);
    }

    // Bland's rule for pivoting.
    // Choose entering index to be lexicographically first with s_j < 0.
    int Simplex::select_entering_variable_Bland(const Eigen::VectorXd & s) const
    {
        const auto M = m_lp.get_num_rows();
        const auto N = m_lp.get_num_vars();
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
        const auto M = m_lp.get_num_rows();
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
        const auto M = m_lp.get_num_rows();
        const auto N = m_lp.get_num_vars();
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
                                         const Basis & basis,
                                         int entering_index) const
    {
        if (m_lp.has_non_trivial_upper_bounds())
        {
            return select_leaving_variable_SUB(x, d, basis, entering_index);
        }
        return select_leaving_variable_Bland(x, d);
    }

    // Simple Upper Bound rule for pivoting.
    int Simplex::select_leaving_variable_SUB(const Eigen::VectorXd & x, const Eigen::VectorXd & d, const Basis & basis,
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
                const auto it = m_lp.var_ubnd.find(basis.get_basic_columns()[i]);
                if (it == m_lp.var_ubnd.end())
                {
                    throw "No upper bound found for basis element: " + std::to_string(basis.get_basic_columns()[i]);
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

        const auto it = m_lp.var_ubnd.find(basis.get_non_basic_columns()[entering_index]);
        double ub_s = it->second;

        if (ub_s < min_theta)
        {
            m_lp.upper_bound_substitution(basis.get_non_basic_columns()[entering_index], ub_s);
            return -2;  // do not perform pivot
        }
        if (need_substitution)
        {
            const auto k = basis.get_basic_columns()[min_i];

            const auto it1 = m_lp.var_ubnd.find(k);
            if (it1 == m_lp.var_ubnd.end())
            {
                throw "No upper bound found for basis element: " + std::to_string(k);
            }

            m_lp.upper_bound_substitution(k, it1->second);
        }

        return min_i;  // perform usual pivot (except if min_i == -1)
    }

    void Simplex::solution_found(Eigen::VectorXd vector_bx, const Eigen::VectorXd & vector_c_B, const Basis & basis)
    {
        double cost = 0.0;
        m_solution.clear();

        const std::vector<Eigen::Index> basis_columns{ basis.get_basic_columns().begin(),
                                                       basis.get_basic_columns().end() };
        for (size_t i = 0u; i < basis_columns.size(); ++i)
        {
            const auto var_id = basis_columns[i];
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

        for (auto var_id : basis.get_non_basic_columns())
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

    int Simplex::select_entering_variable(const Eigen::VectorXd & s) const
    {
        return select_entering_variable_Bland(s);
        // return select_entering_variable_most_neg(s);
    }

    // The (Revised) Simplex Algorithm.
    // arg_basis : on input: initial basis; on result: final basis
    int Simplex::simplex(Basis & basis)
    {
        const auto M = m_lp.get_num_rows();
        const auto N = m_lp.get_num_vars();

        if (basis.get_basic_columns().size() != M)
        {
            throw "Basis size must be equal to M=" + std::to_string(M);
        }

        Eigen::MatrixXd matrix_A_B(M, M);
        Eigen::MatrixXd matrix_A_N(M, N - M);
        Eigen::VectorXd vector_c_B(M);
        Eigen::VectorXd vector_c_N(N - M);

        // Initialize upper-bound substitutions.
        m_lp.ub_substitutions.assign(N, false);

        unsigned long iteration_count = 0L;
        constexpr unsigned long iteration_limit = 1000000000L;
        while (iteration_count < iteration_limit)
        {
            iteration_count++;
            LOG(debug) << "SIMPLEX iteration: " << iteration_count;
            LOG(debug) << "Basis: " << basis.show();

            // Split the matrix A into basis matrix (A_B) and non-basis matrix (A_N).
            // TODO: no need to copy whole vectors, as only 1 element has changed in basis.
            // However note that matrix_A, vector_b and vector_c might have changed
            // due to applying upper-bound substitutions.
            int jB = 0, jN = 0;
            for (auto index : basis.get_basic_columns())
            {
                matrix_A_B.col(jB) = m_lp.matrix_A.col(index);
                vector_c_B[jB] = m_lp.vector_c[index];
                jB++;
            }
            for (auto index : basis.get_non_basic_columns())
            {
                matrix_A_N.col(jN) = m_lp.matrix_A.col(index);
                vector_c_N[jN] = m_lp.vector_c[index];
                jN++;
            }

            //LOG(debug) << "A_B=\n" << matrix_A_B;
            //LOG(debug) << "A_N=\n" << matrix_A_N;
            //LOG(debug) << "c_B=\n" << vector_c_B;
            //LOG(debug) << "c_N=\n" << vector_c_N;

            // Compute x = A_B^{-1} * b and y = (A_B^T)^{-1} * c_B
            Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXd> > lu(matrix_A_B);
            const Eigen::VectorXd x = lu.solve(m_lp.vector_b);
            const Eigen::VectorXd y = lu.transpose().solve(vector_c_B);

            //LOG(debug) << "solved x=\n" << x;
            //LOG(debug) << "solved y=\n" << y;

            // Pricing (reduced costs): s = c_N - (A_N)^T * y
            Eigen::VectorXd s = vector_c_N - matrix_A_N.transpose() * y;

            //LOG(debug) << "reduced costs s=\n" << s;

            // 1) Select the entering variable.
            const int entering_index = select_entering_variable(s);
            if (entering_index == -1)
            {
                solution_found(x, vector_c_B, basis);
                break;
            }

            // Copy the entering column A_q into vector_cy.
            Eigen::VectorXd A_entering = m_lp.matrix_A.col(basis.get_non_basic_columns()[entering_index]);

            // 2) Select the leaving variable.

            // Solve: A_B * d = A(entering_index). Note that matrix_A_B already contains LU factors.
            const Eigen::VectorXd d = lu.solve(A_entering);

            const auto leaving_index = select_leaving_variable(x, d, basis, entering_index);
            if (leaving_index == -1)
            {
                LOG(debug) << "LP IS UNBOUNDED.\n";
                break;
            }

            if (leaving_index != -2)
            {
                // Update basis.
                basis.update(entering_index, leaving_index);
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
