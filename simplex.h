// Copyright (C) 2017-2019 Maciej Drwal
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
// A standard form LP is:
// min c^T x
// s.t.: Ax = b, x >= 0.
//

struct Constraint
{
	char type;
	double rhs;
	map<string, double> name_coeff;     // variable name and its coefficient a_{ij}
	
	Constraint(char _type, double _rhs, map<string, double> _name_coeff) :
		type(_type),
		rhs(_rhs),
		name_coeff(_name_coeff) {}
    
    ~Constraint() = default;
};

struct LinearProgram
{
    int N;
    int M;
	
    map<string, shared_ptr<Constraint> > constraints;
	map<string, double> objective_name_coeff;
    
	LinearProgram() :
		matrix_A(nullptr),
		vector_b(nullptr),
		vector_c(nullptr),
		N(0),
		M(0),
		sense('m'),
		all_inequalities(true),
        objective_label("") {}
        
    ~LinearProgram() = default;
		
	void solve();
    void write(const string & filename) const;
    
    void add_variable(const string & var_name);
    bool has_variable(const string & var_name) const;
    void set_sense(const char s) { sense = s; }
    void set_all_inequalities(bool b) { all_inequalities = b; }
    void set_objective_label(const string & label) { objective_label = label; }
    const string & get_objective_label() const { return objective_label; }
    void add_shift(const string & var_name, double value) { var_shifts[var_name] = value; }
    double get_shift(const string & var_name) { return var_shifts[var_name]; }
    map<string, double> & get_var_shifts() { return var_shifts; }
	
private:
    
    double * matrix_A;
    double * vector_b;
    double * vector_c;
    
	char sense;	// M=maximize / m=minimize
    double objective_value;
	bool all_inequalities;
    string objective_label;
    
	map<string, int> variable_names_to_ids;
    map<string, double> var_shifts;
    map<string, double> solution;
	
	static int * ipiv;

    void presolve();
    void initialize_tableau();
	void print_tableau() const;
    int simplex(set<int> & initial_basis);
	void lineq_solve(double * _matrix_A, double * _vector_bx, double * _vector_cy, bool factorize);
    
};

#endif
