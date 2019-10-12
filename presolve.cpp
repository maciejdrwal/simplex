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

using std::string;
using std::map;
using std::shared_ptr;
using std::cout;
using std::endl;

void Presolve::fix_variable(const string & var_name, double values) {
    throw "Presolve: not implemented.";
}

double Presolve::get_shift(const string & var_name) const
{ 
    return lp->var_shifts[var_name]; 
}

void Presolve::set_shift(const string & var_name, double value)
{ 
    lp->var_shifts[var_name] = value; 
}

// Eliminate all LBs by substituting:
// x' = x - LB
// x' >= 0
void Presolve::eliminate_lbs()
{
    for (map<string, double>::iterator it = lp->var_lbnd.begin(); it != lp->var_lbnd.end(); ++it) 
    {
        string var_name = it->first;
        double lb = it->second;
        auto ub_it = lp->var_ubnd.find(var_name);
        if (ub_it != lp->var_ubnd.end())
        {
            double ub = ub_it->second;

            if (lb > ub) {
                std::cout << "Presolve: problem infeasible, lb>ub for variable: " << var_name << std::endl;
                std::exit(0);
            }
            else if (_isfloatzero(_abs(ub - lb))) {
                // Substitute the variable by its LB;
                fix_variable(var_name, lb);
            }
            else {
                // Update also corresponding UB
                ub_it->second = ub - lb;
            }
        }

        if (_isfloatzero(lb)) continue;

        lp->var_shifts[var_name] = lb;
        it->second = 0.0;

        // Update objective function
        auto obj_fun_it = lp->objective_name_coeff.find(var_name);
        if (obj_fun_it != lp->objective_name_coeff.end()) {
            std::cout << "Presolve: applying shift " << lb << " to obj.fun. variable:" << var_name << std::endl;
            lp->obj_value_shift += (obj_fun_it->second * lb);
        }

        // Update constraints
        for (map<string, shared_ptr<Constraint> >::iterator kt = lp->constraints.begin(); kt != lp->constraints.end(); ++kt) {
            char _type = kt->second->type;
            double _rhs = kt->second->rhs;
            auto term_it = kt->second->name_coeff.find(var_name);
            if (term_it != kt->second->name_coeff.end()) {
                std::cout << "Presolve: applying shift " << lb << " to constraint " << kt->first 
                    << " variable " << var_name << std::endl;
                kt->second->rhs -= (lb * term_it->second);
            }
        }
    }
}

void Presolve::apply_reductions()
{
    for (auto it = lp->constraints.begin(); it != lp->constraints.end(); ++it) {
        double U = 0.0;
        double L = 0.0;
        for (auto nc_it = it->second->name_coeff.begin(); nc_it != it->second->name_coeff.end(); ++nc_it) {
            double a = nc_it->second;
            double ub = lp->var_ubnd[nc_it->first];
            double lb = lp->var_lbnd[nc_it->first];
            if (nc_it->second > 0.0) {
                U += a * ub;
                L += a * lb;
            }
            else if (nc_it->second < 0.0) {
                U += a * lb;
                L += a * ub;
            }
        }
        cout << "apply_red: constr:" << it->first << " L= " << L << " U=" << U << endl;
    }
}

void Presolve::run() 
{
    if (reductions_enabled) apply_reductions();
    eliminate_lbs();
}

