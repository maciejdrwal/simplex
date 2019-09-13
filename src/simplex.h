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
	map<string, double> name_coeff;
	
	Constraint(char _type, double _rhs, map<string, double> _name_coeff) :
		type(_type),
		rhs(_rhs),
		name_coeff(_name_coeff) {}
};

struct LinearProgram
{
    double * matrix_A;
    double * vector_b;
    double * vector_c;
    
    int N;
    int M;
	
	char sense;	// M=maximize / m=minimize
    
    string objective_label;
    double objective_value;     // Placeholder for optimal value.
	map<string, int> variable_names;        
	map<string, double> objective_name_coeff;
	map<string, double> var_shifts;
    map<string, double> solution;
	bool all_inequalities;	
	
	map<string, Constraint*> constraints;
    // map<string, shared_ptr<Constraint> > constraints;
    
	LinearProgram() :
		matrix_A(NULL),
		vector_b(NULL),
		vector_c(NULL),
		N(0),
		M(0),
		sense('m'),
		all_inequalities(true) {}
	
	~LinearProgram() {
		// Free constraints memory.
		for (map<string, Constraint*>::iterator it = constraints.begin(); it != constraints.end(); ++it) {
			delete it->second;
		}
		constraints.clear();
	}
	
	void solve();
    void write(const string & filename);
	
private:
	
	static int * ipiv;

    void initialize_tableau();
	void print_tableau();
    int simplex(set<int> & initial_basis);
	void lineq_solve(double * _matrix_A, double * _vector_bx, double * _vector_cy, bool factorize);
    
};

#endif
