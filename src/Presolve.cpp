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

#include <map>
#include <set>

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
        std::set<std::string> redundant;
        for (auto & [label, constraint] : m_lp.constraints)
        {
            double U = 0.0;
            double L = 0.0;
            for (const auto & [var_name, a] : constraint.name_to_coeff)
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
            LOG(debug) << "Presolve: applying reduction to constraint:" << label << " L=" << L << " U=" << U;
            if ((constraint.type == '<' && U <= constraint.rhs) || (constraint.type == '>' && L >= constraint.rhs))
            {
                LOG(debug) << "  constraint is redundant";
                redundant.insert(label);
            }
            if (((constraint.type == '<' || constraint.type == '=') && L > constraint.rhs) ||
                ((constraint.type == '>' || constraint.type == '=') && U < constraint.rhs))
            {
                LOG(debug) << "  no feasible solution";
                std::exit(0);
            }
        }

        for (const auto & label : redundant)
        {
            m_lp.constraints.erase(label);
        }
    }

    void Presolve::run()
    {
        if (m_reductions_enabled)
        {
            apply_reductions();
        }
        eliminate_lbs();
    }
}