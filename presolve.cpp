#include <map>

#include "presolve.h"

#define _abs(x) ((x) < 0.0 ? -(x) : (x))
#define _isfloatzero(x) (_abs(x) < 1e-9 ? true : false)

void apply_shifts(LinearProgram * lp) {
    // printf("presolve: variables -> indices:\n");
    // for (map<string, int>::iterator it = lp->variable_names.begin(); it != lp->variable_names.end(); ++it) {
    //     printf("%s -> %i\n", it->first.c_str(), it->second);
    // }
    
    for (map<string, double>::iterator it = lp->objective_name_coeff.begin(); it != lp->objective_name_coeff.end(); ++it) {
        // Applying shifts.
        for (map<string, double>::iterator kt = lp->get_var_shifts().begin(); kt != lp->get_var_shifts().end(); ++kt) {
            if ((kt->first == it->first) && (_isfloatzero(kt->second))) {
                printf("APPLYING ZERO SHIFT: %s\n", it->first.c_str());
                it->second = -it->second;
            }
        }
        
        printf("%f * %s + ", it->second, it->first.c_str());
    }
    printf("\n\nConstraints:\n");
    for (map<string, shared_ptr<Constraint> >::iterator it = lp->constraints.begin(); it != lp->constraints.end(); ++it) {
        char _type = it->second->type;
        double _rhs = it->second->rhs;
        if ((_type == '<' && _rhs < 0.0) ||
            (_type == '>' && _rhs > 0.0) ||
            (_type == '=')) lp->set_all_inequalities(false);

        printf("%s : ", it->first.c_str());
        for (map<string, double>::iterator jt = it->second->name_coeff.begin(); jt != it->second->name_coeff.end(); ++jt) {
            printf("%f %s + ", jt->second, jt->first.c_str());
            
            // Applying shifts.     
            for (map<string, double>::iterator kt = lp->get_var_shifts().begin(); kt != lp->get_var_shifts().end(); ++kt) {
                if (kt->first == jt->first) {
                    if (_isfloatzero(kt->second)) {
                        printf("APPLYING ZERO SHIFT: %s -> %f\n", kt->first.c_str(), kt->second);
                        jt->second = -jt->second;
                    }
                    else {
                        printf("APPLYING POS. SHIFT: %s -> %f\n", kt->first.c_str(), kt->second);
                        it->second->rhs -= kt->second * jt->second;
                    }
                }       
            }
        }
        printf("%c %f\n", it->second->type, it->second->rhs);
    }
}
