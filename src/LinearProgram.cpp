// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#include "LinearProgram.h"

#include "Logger.h"
#include "Utils.h"

namespace
{
    const std::string ARTIFICIAL = "_ARTIFICIAL_";
}

namespace simplex
{
    void Constraint::add_term(std::string_view name, double a)
    {
        const auto it = name_to_coeff.find(name.data());
        if (it == name_to_coeff.end())
        {
            name_to_coeff[name.data()] = a;  // add new term
        }
        else
        {
            it->second += a;  // add to existing term
        }
    }

    bool Constraint::has_variable(std::string_view name) const
    { return name_to_coeff.count(name.data()) > 0; }

    std::optional<double> Constraint::get_coefficient(std::string_view name) const
    {
        const auto it = name_to_coeff.find(name.data());
        return it == name_to_coeff.end() ? std::nullopt : std::optional(it->second);
    }

    void Constraint::remove_term(std::string_view name)
    {
        const auto it = name_to_coeff.find(name.data());
        if (it != name_to_coeff.end())
        {
            name_to_coeff.erase(it);
        }
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
                    constraint.add_term(var_name);
                }
                if (constraint.type == '>')
                {
                    constraint.add_term(var_name, -1.0);
                    all_inequalities = false;
                }
                constraint.type = '=';
                const auto var_id = add_variable(var_name);

                // If all constraints were inequalities then use slacks for initial basis.
                if (all_inequalities)
                {
                    basis.push_back(var_id);
                }
                LOG(debug) << "added slack variable: " << var_name << " (" << var_id << ")";
            }
            else
            {
                all_inequalities = false;
            }
        }
        LOG(debug) << "added " << num_of_slack_vars << " slack variables, all_inequalities=" << all_inequalities;
        if (!all_inequalities)
        {
            basis.clear();
        }
        return all_inequalities;
    }

    int LinearProgram::add_artificial_variables_for_first_phase(Basis & basis)
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
            basis.push_back(var_id);
            LOG(debug) << "added artificial variable: " << var_name;
        }
        return art_var_id;
    }

    std::string LinearProgram::get_artificial_variable(int i)
    {
		return ARTIFICIAL + utils::to_str<int>(i);
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