// Copyright (C) 2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#include <iostream>
#include <map>

#include "presolve.h"
#include "simplex.h"
#include "utils.h"
#include "logger.h"

namespace simplex
{
    void Presolve::fix_variable(const std::string & var_name, double values)
    {
        throw "Presolve: not implemented.";
    }

    double Presolve::get_shift(const std::string & var_name) const
    {
        return m_lp.var_shifts[var_name];
    }

    void Presolve::set_shift(const std::string & var_name, double value)
    {
        m_lp.var_shifts[var_name] = value;
    }

    // Eliminate all LBs by substituting:
    // x' = x - LB
    // x' >= 0
    void Presolve::eliminate_lbs()
    {
        for (auto & [var_name, lb_ref] : m_lp.var_lbnd)
        {
            double lb = lb_ref;
            if (utils::isfloatzero(lb)) continue;

            auto ub_it = m_lp.var_ubnd.find(var_name);
            if (ub_it != m_lp.var_ubnd.end())
            {
                double ub = ub_it->second;

                if (lb > ub)
                {
                    LOG(debug) << "Presolve: problem infeasible, lb>ub for variable: " << var_name;
                    std::exit(0);
                }
                else if (utils::isfloatzero(utils::abs(ub - lb)))
                {
                    // Substitute the variable by its LB;
                    fix_variable(var_name, lb);
                }
                else
                {
                    // Update also corresponding UB
                    ub_it->second = ub - lb;
                }
            }

            m_lp.var_shifts[var_name] = lb;
            lb_ref = 0.0;

            // Update objective function
            auto obj_fun_it = m_lp.objective_name_coeff.find(var_name);
            if (obj_fun_it != m_lp.objective_name_coeff.end())
            {
                LOG(debug) << "Presolve: applying shift " << lb << " to obj.fun. variable:" << var_name;
                m_lp.obj_value_shift += (obj_fun_it->second * lb);
            }

            // Update constraints
            for (auto & [constr_name, constraint] : m_lp.constraints)
            {
                const auto a = constraint.get_coefficient(var_name);
                if (a.has_value())
                {
                    LOG(debug) << "Presolve: applying shift " << lb << " to constraint " << constr_name << " variable " << var_name;
                    constraint.rhs -= (lb * a.value());
                }
            }
        }
    }

    void Presolve::apply_reductions()
    {
        for (const auto & constraint : m_lp.constraints)
        {
            const auto & data = constraint.second;
            double U = 0.0;
            double L = 0.0;
            for (const auto & [var_name, a] : data.name_to_coeff)
            {
                // coefficient a_ij * x_i
                double ub = m_lp.var_ubnd[var_name];
                double lb = m_lp.var_lbnd[var_name];
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
            LOG(debug) << "apply_reductions: constr:" << constraint.first << " L=" << L << " U=" << U;
            if ((data.type == '<' && U <= data.rhs) ||
                (data.type == '>' && L >= data.rhs))
            {
                LOG(debug) << "  constraint is redundant";
            }
            if (((data.type == '<' || data.type == '=') && L > data.rhs) ||
                ((data.type == '>' || data.type == '=') && U < data.rhs))
            {
                LOG(debug) << "  no feasible sol.";
            }
        }
    }

    void Presolve::run()
    {
        if (b_reductions_enabled)
        {
            apply_reductions();
        }
        eliminate_lbs();
    }
}
