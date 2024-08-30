// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#include "Presolve.h"
#include "Simplex.h"
#include "Utils.h"
#include "Logger.h"

#include <iostream>
#include <map>

namespace simplex
{
    void Presolve::fix_variable(Eigen::Index var_id, double values)
    {
        throw "Presolve: not implemented.";
    }

    // Eliminate all LBs by substituting:
    // x' = x - LB
    // x' >= 0
    void Presolve::eliminate_lbs()
    {
        for (auto & [var_id, lb_ref] : m_lp.var_lbnd)
        {
            if (utils::is_float_zero(lb_ref))
            {
                continue;
            }

            double lb = lb_ref;
            auto ub_it = m_lp.var_ubnd.find(var_id);
            if (ub_it != m_lp.var_ubnd.end())
            {
                double ub = ub_it->second;

                if (lb > ub)
                {
                    LOG(debug) << "Presolve: problem infeasible, lb>ub for variable: " << var_id;
                    std::exit(0);
                }
                else if (utils::is_float_zero(utils::abs(ub - lb)))
                {
                    // Substitute the variable by its LB
                    fix_variable(var_id, lb);
                }
                else
                {
                    // Update also corresponding UB
                    ub_it->second = ub - lb;
                }
            }

            m_lp.var_shifts[var_id] = lb;
            lb_ref = 0.0;

            const auto & var_name = m_lp.variable_id_to_name[var_id];

            // Update objective function
            auto obj_fun_it = m_lp.objective_coeff.find(var_id);
            if (obj_fun_it != m_lp.objective_coeff.end())
            {
                LOG(debug) << "Presolve: applying shift " << lb << " to obj.fun. variable: " << var_name;
                m_lp.obj_value_shift += (obj_fun_it->second * lb);
            }

            // Update constraints
            for (auto & [constr_name, constraint] : m_lp.constraints)
            {
                const auto a = constraint.get_coefficient(var_name);
                if (a.has_value())
                {
                    LOG(debug) << "Presolve: applying shift " << lb << " to constraint " << constr_name << " variable: " << var_name;
                    constraint.rhs -= (lb * a.value());
                }
            }
        }
    }

    void Presolve::apply_reductions() const
    {
        for (const auto & constraint : m_lp.constraints)
        {
            const auto & data = constraint.second;
            double U = 0.0;
            double L = 0.0;
            for (const auto & [var_name, a] : data.name_to_coeff)
            {
                // coefficient a_ij * x_i
                const auto it = m_lp.variable_name_to_id.find(var_name);
                const auto var_id = it->second;
                double ub = m_lp.var_ubnd[var_id];
                double lb = m_lp.var_lbnd[var_id];
                if (a > 0.0)
                {
                    U += a * ub;
                    L += a * lb;
                }
                else if (a < 0.0)
                {
                    U += a * lb;
                    L += a * ub;
                }
            }
            LOG(debug) << "Presolve: applying reduction to constraint:" << constraint.first << " L=" << L << " U=" << U;
            if ((data.type == '<' && U <= data.rhs) ||
                (data.type == '>' && L >= data.rhs))
            {
                LOG(debug) << "  constraint is redundant";
            }
            if (((data.type == '<' || data.type == '=') && L > data.rhs) ||
                ((data.type == '>' || data.type == '=') && U < data.rhs))
            {
                LOG(debug) << "  no feasible solution";
                std::exit(0);
            }
        }
    }

    void Presolve::run()
    {
        if (m_reductions_enabled)
        {
            apply_reductions();
        }
        eliminate_lbs();

        // Add constraint for upper bound
        //     for (auto & [var_id, ub] : m_lp.var_ubnd)
        //     {
        //         const auto it = m_lp.variable_id_to_name.find(var_id);
        //         if (it == m_lp.variable_id_to_name.end())
        //         {
        //	throw "Presolve: variable not found.";
        //}
        //         Constraint constraint('<', ub);
        //         constraint.add_term(it->second);
        //         m_lp.constraints.emplace("UB_" + it->second, std::move(constraint));
        //     }
    }
}