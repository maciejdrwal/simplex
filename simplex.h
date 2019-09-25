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

using namespace std;

//
// A standard form LP:
// min c^T x
// s.t.: Ax = b, 0 <= x <= UB
//

struct Constraint
{
    char type;
    double rhs;
    map<string, double> name_coeff;     // variable name and its coefficient a_{ij}
    
    Constraint(char _type = '<', double _rhs = 0.0) :
        type(_type), rhs(_rhs) {}
    
    ~Constraint() = default;
};

struct LinearProgram
{
    int N;
    int M;
    
    map<string, shared_ptr<Constraint> > constraints;
    map<string, double> objective_name_coeff;
    map<string, double> var_lbnd;   // LBs are set to 0 by variable substitutions in presovle
    map<string, double> var_ubnd;
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
        _var_coeff_cache(1.0),
        _sign_cache(1.0) 
        {}
        
    ~LinearProgram() = default;
        
    void solve();
    void write(const string & filename) const;
    
    void add_variable(const string & var_name);
    bool has_variable(const string & var_name) const;
    void set_sense(const char s) { sense = s; }
    void set_all_inequalities(bool b) { all_inequalities = b; }
    void set_objective_label(const string & label) { objective_label = label; }
    const string & get_objective_label() const { return objective_label; }
    void set_obj_value_shift(double x) { obj_value_shift = x; }
    void set_var_coeff_cache(double x) { _var_coeff_cache = x; }
    double get_var_coeff_cache() const { return _var_coeff_cache; }
    void set_sign_cache(double x) { _sign_cache = x; }
    double get_sign_cache() const { return _sign_cache; }
    void set_label_cache(const string & s) { _label_cache = s; }
    const string & get_label_cache() const { return _label_cache; }

    void set_shift(const string & var_name, double value) { var_shifts[var_name] = value; }
    double get_shift(const string & var_name) { return var_shifts[var_name]; }
    map<string, double> & get_var_shifts() { return var_shifts; }
    
private:
    
    double * matrix_A;
    double * vector_b;
    double * vector_c;
    
    char sense; // M=maximize / m=minimize
    double objective_value;
    string objective_label;
    
    bool all_inequalities;

    // variable cache (used by parser)
    double _var_coeff_cache;
    double _sign_cache;
    string _label_cache;

    vector<int> basis;
    vector<int> non_basis;
    vector<bool> ub_substitutions;

    map<string, int> variable_names_to_ids;
    map<int, string> variable_ids_to_names;
    map<string, double> var_shifts;
    map<string, double> solution;
    
    static int * ipiv;

    void initialize_tableau();
    void print_tableau() const;
    int simplex(set<int> & initial_basis);
    void lineq_solve(double * _matrix_A, double * _vector_bx, double * _vector_cy, bool factorize);
    int select_entering_variable_Bland(double * vector_c_N);
    int select_entering_variable_most_neg(double * vector_c_N);
    int select_leaving_variable_Bland(double * vector_bx, double * vector_cy, double * matrix_A_B);
    int select_leaving_variable_GUB(double * vector_bx, double * vector_cy, double * matrix_A_B, int entering_index);
    void upper_bound_substitution(int var_id, double ub);
    void solution_found(double * vector_bx, double * vector_cy, double * vector_c_B);
};

#endif
