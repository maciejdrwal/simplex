// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#include "LinearProgram.h"

#include "Logger.h"
#include "Utils.h"

namespace simplex
{
    Constraint & Constraint::add_term(std::string_view name, double a)
    {
        if (utils::is_float_zero(a))
        {
            return *this;
        }
        const auto it = name_to_coeff.find(name.data());
        if (it == name_to_coeff.end())
        {
            name_to_coeff[name.data()] = a;  // add new term
        }
        else
        {
            it->second += a;  // add to existing term
        }
        return *this;
    }

    bool Constraint::has_variable(std::string_view name) const
    { return name_to_coeff.count(name.data()) > 0; }

    std::optional<double> Constraint::get_coefficient(std::string_view name) const
    {
        const auto it = name_to_coeff.find(name.data());
        return it == name_to_coeff.end() ? std::nullopt : std::optional(it->second);
    }

    Constraint & Constraint::remove_term(std::string_view name)
    {
        const auto it = name_to_coeff.find(name.data());
        if (it != name_to_coeff.end())
        {
            name_to_coeff.erase(it);
        }
        return *this;
    }

    void Constraint::negate_sides()
    {
        rhs = -rhs;
        for (auto & [_, coeff] : name_to_coeff)
        {
            coeff = -coeff;
        }
        type = (type == '<') ? '>' : ((type == '>') ? '<' : type);
    }

    // This should be called just before simplex(), after presolve and all
    // other problem transformations.
    void LinearProgram::initialize_tableau()
    {
        const auto M = get_num_rows();
        const auto N = get_num_vars();

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

        // Prepare the original objective function.
        for (const auto [i, coeff] : objective_coeff)
        {
            vector_c[i] = (sense == 'm') ? coeff : -coeff;
        }

        rescale_matrix();
    }

    void LinearProgram::print_tableau() const
    {
        // Note that matrix_A is stored in a column-major order.
        std::stringstream ss;
        ss << "A=\n" << matrix_A << '\n';
        ss << "b=\n" << vector_b << '\n';
        ss << "c=\n" << vector_c << '\n';
        LOG(debug) << ss.str();
    }


    size_t LinearProgram::add_slack_variables_for_inequality_constraints()
    {
        size_t num_of_slack_vars = 0;
        bool all_inequalities = true;
        for (auto & [label, constraint] : constraints)
        {
            if (constraint.type != '=')
            {
                const std::string var_name("__SLACK" + std::to_string(num_of_slack_vars++));
                if (constraint.type == '<')
                {
                    constraint.add_term(var_name);
                }
                if (constraint.type == '>')
                {
                    constraint.negate_sides();
                    constraint.add_term(var_name);
                }
                constraint.type = '=';
                const auto var_id = add_variable(var_name);
                m_slacks[label] = var_id;

                LOG(debug) << "added slack variable: " << var_name << " (" << var_id << ")";
            }
            else
            {
                all_inequalities = false;
            }
        }
        LOG(debug) << "added " << num_of_slack_vars << " slack variables, all_inequalities=" << all_inequalities;
        return num_of_slack_vars;
    }

    void LinearProgram::set_lower_bound(Eigen::Index var_id, double low_value)
    {
        var_lbnd[var_id] = low_value;
    }

    void LinearProgram::set_upper_bound(Eigen::Index var_id, double low_value)
    {
        var_ubnd[var_id] = low_value;
        m_has_UBS = true;
    }

    void LinearProgram::upper_bound_substitution(Eigen::Index var_id, double ub)
    {
        // Modify the objective function: c_s x_s -> -c_s x_s
        vector_c[var_id] = -vector_c[var_id];

        // Modify constraints: a_{js} x_s -> -a_{js} x_s, b_j -> b_j - a_{js}*u_s
        const auto M = constraints.size();
        for (auto i = 0u; i < M; i++)
        {
            const auto a = matrix_A(i, var_id);
            vector_b[i] = vector_b[i] - a * ub;
            matrix_A(i, var_id) = -a;
        }

        ub_substitutions[var_id] = !ub_substitutions[var_id];
    }

    void LinearProgram::rescale_matrix()
    {
        const auto N = matrix_A.rows();
        m_row_scaling_factors.resize(N);
        for (auto j = 0u; j < N; j++)
        {
            const auto max_a = matrix_A.row(j).cwiseAbs().maxCoeff();
            m_row_scaling_factors[j] = max_a;
            if (max_a > 1.0)
            {
                matrix_A.row(j) /= max_a;
                vector_b[j] /= max_a;
            }
        }

        const auto M = matrix_A.cols();
        m_column_scaling_factors.resize(M);
        for (auto j = 0u; j < M; j++)
        {
            const auto max_a = matrix_A.col(j).cwiseAbs().maxCoeff();
            m_column_scaling_factors[j] = max_a;
            matrix_A.col(j) /= max_a;
            vector_c[j] /= max_a;
        }
    }

    void LinearProgram::add_constraint(const std::string & name, Constraint && constraint)
    {
        constraints.emplace(name, std::move(constraint));
    }

    const Constraint & LinearProgram::get_constraint(const std::string & name) const
    {
        const auto it = constraints.find(name);
        if (it != constraints.end())
        {
            return it->second;
        }
        throw "Constraint not found: " + name;
    }

    // Note: each variable has ID that corresponds to that variable's position
    // in the computation array vector_c, and is also used, e.g., to
    // indicate basic/non-basic variables. 
    ///@param var_name - name of the new variable
    ///@param coeff - the coefficient of the variable in the objective function
    ///@return the index of the new variable
    Eigen::Index LinearProgram::add_variable(const std::string & var_name, double coeff)
    {
        const auto var_id = next_variable_id++;
        variable_name_to_id[var_name] = var_id;
        variable_id_to_name[var_id] = var_name;
        objective_coeff[var_id] = coeff;
        var_lbnd[var_id] = 0.0;  // set default bounds
        var_ubnd[var_id] = std::numeric_limits<double>::max();
        return var_id;
    }

    void LinearProgram::remove_variable(Eigen::Index var_id)
    {
        const auto it = variable_id_to_name.find(var_id);
        if (it == variable_id_to_name.end())
        {
            return;
        }

        objective_coeff.erase(var_id);
        var_lbnd.erase(var_id);
        var_ubnd.erase(var_id);
        variable_name_to_id.erase(it->second);
        variable_id_to_name.erase(var_id);
    }

    bool LinearProgram::has_variable(const std::string & var_name) const
    {
        return variable_name_to_id.count(var_name) > 0;
    }

    Eigen::Index LinearProgram::get_variable_id(const std::string & var_name) const
    {
        const auto it = variable_name_to_id.find(var_name);
        if (it == variable_name_to_id.end())
        {
            throw "Variable not found: " + var_name;
        }
        return it->second;
    }
}