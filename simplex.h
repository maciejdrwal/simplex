// Copyright (C) 2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#include "presolve.h"

struct Constraint
{
    char type;
    double rhs;
    std::map<std::string, double> name_coeff;     // variable name and its coefficient a_{ij}
    
    Constraint(char _type = '<', double _rhs = 0.0) :
        type(_type), rhs(_rhs) {}
    
    ~Constraint() = default;
};

//
// A standard form LP:
// min c^T x
// s.t.: Ax = b, 0 <= x <= UB
//
struct LinearProgram
{
    friend class Presolve;

    int N;
    int M;
    
    std::map<std::string, std::shared_ptr<Constraint> > constraints;
    std::map<std::string, double> objective_name_coeff;   // maps variable names x_j to their c_j
    std::map<std::string, double> var_lbnd;   // LBs are set to 0 by variable substitutions in presovle
    std::map<std::string, double> var_ubnd;   // UBs are modified by presolve
    double obj_value_shift;         // constant term in objective function
    
    LinearProgram() :
        matrix_A(nullptr),
        vector_b(nullptr),
        vector_c(nullptr),
        N(0),
        M(0),
        sense('m'),
        all_inequalities(true),
        objective_label(""),
        objective_value(0.0),
        obj_value_shift(0.0),
        presolve(this, true)
        {}
        
    ~LinearProgram() = default;
        
    void solve();
    void write(const std::string & filename) const;
    
    void add_variable(const std::string & var_name);
    void remove_variable(const std::string & var_name);
    bool has_variable(const std::string & var_name) const;
    void set_sense(const char s) { sense = s; }
    void set_all_inequalities(bool b) { all_inequalities = b; }
    void set_objective_label(const std::string & label) { objective_label = label; }
    const std::string & get_objective_label() const { return objective_label; }
        
private:
    
    std::unique_ptr<double[]> matrix_A;
    std::unique_ptr<double[]> vector_b;
    std::unique_ptr<double[]> vector_c;
    
    char sense; // M=maximize / m=minimize
    double objective_value;
    std::string objective_label;
    
    bool all_inequalities;

    std::vector<int> basis;
    std::vector<int> non_basis;
    std::vector<bool> ub_substitutions;

    std::map<std::string, int> variable_name_to_id;
    std::map<int, std::string> variable_id_to_name;
    std::map<std::string, double> var_shifts;
    std::map<std::string, double> solution;

    Presolve presolve;
    
    static int * ipiv;

    void initialize_tableau();
    void print_tableau() const;
    int simplex(std::set<std::string> & initial_basis);
    void lineq_solve(double * _matrix_A, double * _vector_bx, double * _vector_cy, bool factorize);
    int select_entering_variable_Bland(double * vector_c_N);
    int select_entering_variable_most_neg(double * vector_c_N);
    int select_leaving_variable_Bland(double * vector_bx, double * vector_cy, double * matrix_A_B);
    int select_leaving_variable_SUB(double * vector_bx, double * vector_cy, double * matrix_A_B, int entering_index);
    void upper_bound_substitution(int var_id, double ub);
    void solution_found(double * vector_bx, double * vector_cy, double * vector_c_B);
};

#endif
