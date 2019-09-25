// Copyright (C) 2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#include <iostream>
#include <map>

#include "presolve.h"
#include "utils.h"

void fix_variable(const string & var_name, double values) {
    throw "Presolve: not implemented.";
}

void eliminate_lbs(LinearProgram * lp)
{
    // Eliminate LBs by substituting:
    // z = x - LB
    // z >= 0
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

        lp->set_shift(var_name, lb);
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

void presolve(LinearProgram * lp) 
{
    eliminate_lbs(lp);
}

// void apply_shifts(LinearProgram * lp) {
//     // printf("presolve: variables -> indices:\n");
//     // for (map<string, int>::iterator it = lp->variable_names.begin(); it != lp->variable_names.end(); ++it) {
//     //     printf("%s -> %i\n", it->first.c_str(), it->second);
//     // }
    
//     for (map<string, double>::iterator it = lp->objective_name_coeff.begin(); it != lp->objective_name_coeff.end(); ++it) {
//         // Applying shifts.
//         for (map<string, double>::iterator kt = lp->get_var_shifts().begin(); kt != lp->get_var_shifts().end(); ++kt) {
//             if ((kt->first == it->first) && (_isfloatzero(kt->second))) {
//                 printf("APPLYING ZERO SHIFT: %s\n", it->first.c_str());
//                 it->second = -it->second;
//             }
//         }
        
//         printf("%f * %s + ", it->second, it->first.c_str());
//     }
//     printf("\n\nConstraints:\n");
//     for (map<string, shared_ptr<Constraint> >::iterator it = lp->constraints.begin(); it != lp->constraints.end(); ++it) {
//         char _type = it->second->type;
//         double _rhs = it->second->rhs;
//         if ((_type == '<' && _rhs < 0.0) ||
//             (_type == '>' && _rhs > 0.0) ||
//             (_type == '=')) lp->set_all_inequalities(false);

//         printf("%s : ", it->first.c_str());
//         for (map<string, double>::iterator jt = it->second->name_coeff.begin(); jt != it->second->name_coeff.end(); ++jt) {
//             printf("%f %s + ", jt->second, jt->first.c_str());
            
//             // Applying shifts.     
//             for (map<string, double>::iterator kt = lp->get_var_shifts().begin(); kt != lp->get_var_shifts().end(); ++kt) {
//                 if (kt->first == jt->first) {
//                     if (_isfloatzero(kt->second)) {
//                         printf("APPLYING ZERO SHIFT: %s -> %f\n", kt->first.c_str(), kt->second);
//                         jt->second = -jt->second;
//                     }
//                     else {
//                         printf("APPLYING POS. SHIFT: %s -> %f\n", kt->first.c_str(), kt->second);
//                         it->second->rhs -= kt->second * jt->second;
//                     }
//                 }       
//             }
//         }
//         printf("%c %f\n", it->second->type, it->second->rhs);
//     }
// }
