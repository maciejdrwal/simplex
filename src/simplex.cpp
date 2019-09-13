// Copyright (C) 2017-2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#include <iterator>
#include <limits>
#include <chrono>
#include <fstream>

#include <Accelerate/Accelerate.h>

#include "simplex.h"
#include "utils.h"

int * LinearProgram::ipiv = NULL;

chrono::time_point<chrono::steady_clock> t_start, t_end;
map<int, string> variable_names_rev;

// Print internal data structures on screen (for debugging).
void LinearProgram::print_tableau()
{
	// Note that matrix_A is stored in a column-major order.
	printf("N=%d, M=%d\n", N, M);
	printf("A=\n");
	for (int i = 0; i < M; i++) {
	    for (int j = 0; j < N; j++) {
	        printf("%f\t", matrix_A[j * M + i]);
	    }
	    printf("\n");
	}
	printf("b=\n");
	for (int i = 0; i < M; i++) {
	   printf("%f\t", vector_b[i]);
	}
	printf("\nc=\n");
	for (int i = 0; i < N; i++) {
	   printf("%f\t", vector_c[i]);
	}
	printf("\n");
}

// Copies the data read by parser into internal data structures.
void LinearProgram::initialize_tableau()
{
	int j = 0;
	for (map<string, Constraint*>::iterator it = constraints.begin(); it != constraints.end(); ++it) {
		for (map<string, double>::iterator jt = it->second->name_coeff.begin(); jt != it->second->name_coeff.end(); ++jt) {
			int i = variable_names[jt->first] * M + j;
			matrix_A[i] = jt->second;
		}
		vector_b[j++] = it->second->rhs;
	}
	//print_tableau();
}

inline void measure_time_start()
{
	t_start = chrono::steady_clock::now();
}

inline void print_elapsed_time()
{
	t_end = chrono::steady_clock::now();
	printf("Elapsed time: %lld s (%lld ms, %lld us, %lld ns)\n", 
	chrono::duration_cast<chrono::seconds>(t_end - t_start).count(), 
	chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count(),
	chrono::duration_cast<chrono::microseconds>(t_end - t_start).count(),
	chrono::duration_cast<chrono::nanoseconds>(t_end - t_start).count());
}

// The main procedure of the solver.
// Initializes internal data structures and runs the solution algorithm.
// TODO: refactor
void LinearProgram::solve()
{
	measure_time_start();
			
	M = constraints.size();
	N = variable_names.size();
	
	printf("Allocating memory...\n");
	vector_b = new double[M];
	vector_c = new double[N + 2 * M];	// add reserve for slack/artificial variables
	matrix_A = new double[M * (N + 2 * M)];
	
	memset(vector_b, 0, M * sizeof(double));
	memset(vector_c, 0, (N + 2 * M) * sizeof(double));
	memset(matrix_A, 0, (M * (N + 2 * M)) * sizeof(double));
				
	set<int> basis;		// Initial basis.
    
    // Create reverse mapping for variable ids.    
    variable_names_rev.clear();
    for (map<string, int>::iterator it = variable_names.begin(); it != variable_names.end(); ++it) {
        variable_names_rev[it->second] = it->first;
    }
	
	// Add slack variables to inequality constraints.
	int num_of_slack_vars = 0;
	for (map<string, Constraint*>::iterator it = constraints.begin(); it != constraints.end(); ++it) {
		if (it->second->type != '=') {
			string var_name = string("__SLACK") + tostr<int>(num_of_slack_vars++);
			if (it->second->type == '<') {
				it->second->name_coeff[var_name] = 1.0;
			}
			if (it->second->type == '>') {
				it->second->name_coeff[var_name] = -1.0;
			}
			it->second->type = '=';
			variable_names[var_name] = variable_names.size();
			
			// If all constraints are inequalities then use slacks for initial basis.
			if (all_inequalities) basis.insert(variable_names[var_name]);
			printf("added slack variable:%s\n", var_name.c_str());
		}
	}
	int num_of_primal_vars = objective_name_coeff.size();
	
	N = variable_names.size();	// this might have been modified by adding new variables
	
	if (!all_inequalities) {
		// Construct artificial variables for Phase I.
		int i = 0;
		printf("Not all constrains are inequalities: adding artificial variables.\n");
		for (map<string, Constraint*>::iterator it = constraints.begin(); it != constraints.end(); ++it) {
			if (it->second->rhs < 0.0) {
				it->second->rhs = -it->second->rhs;
				for (map<string, double>::iterator jt = it->second->name_coeff.begin(); jt != it->second->name_coeff.end(); ++jt) {
					jt->second = -jt->second;
				}
				// TODO: at this point there should be no inequality constraints
				//if (it->second->type == '<') it->second->type = '>';
				//else if (it->second->type == '>') it->second->type = '<';
				if (it->second->type == '<' || it->second->type == '>') {
					printf("FATAL ERROR WHILE ADDING ART.VARS.\n");
					exit(1);
				}
			}
			string var_name = string("__ARTIFICIAL") + tostr<int>(i++);
			it->second->name_coeff[var_name] = 1.0;
			variable_names[var_name] = variable_names.size();
			basis.insert(variable_names[var_name]);
			printf("added artificial variable:%s\n", var_name.c_str());
		}
		for (int i = 0; i < M; i++) vector_c[num_of_primal_vars + num_of_slack_vars + i] = 1.0;
		N = variable_names.size();	// this has changed after adding artificial variables
		
		// Solve Phase I LP.
		printf("*** PHASE I ***\n");
		initialize_tableau();	
		simplex(basis);
		
		// Remove artificial variables.
		i = 0;
		for (map<string, Constraint*>::iterator it = constraints.begin(); it != constraints.end(); ++it) {
			string var_name = string("__ARTIFICIAL") + tostr<int>(i++);
			
			// Check if some artificial variable remained in the basis.
			if (basis.count(variable_names[var_name]) > 0) {
				// TODO: handle this case
				
				basis.erase(variable_names[var_name]);
				for (map<string, int>::iterator jt = variable_names.begin(); jt != variable_names.end(); ++jt) {
					if (basis.count(jt->second) == 0 && (strncmp(jt->first.c_str(), "__ARTIFICIAL", 12) != 0)) {
						basis.insert(jt->second);
						printf("replaced variable %d (%s) by %d (%s) in the basis\n",  variable_names[var_name], var_name.c_str(), jt->second, jt->first.c_str());
						break;
					}
				}
				
				//printf("Artificial variable %s in the basis. Aborting.\n", var_name.c_str());
				//exit(1);
			}
			
			it->second->name_coeff.erase(var_name);
			variable_names.erase(var_name);
			printf("removed artificial variable:%s\n", var_name.c_str());
		}
		N = variable_names.size();   // this has changed after removing artificial variables
	}
	else {
		printf("Using slacks for initial basis, skipping Phase I.\n");
		initialize_tableau();
	}
	
	if (objective_value > 0.0) {
		printf("LP is infeasible. Aborting.\n");
	}
	else {
		// Prepare the original objective function.
	    for (map<string, double>::iterator it = objective_name_coeff.begin(); it != objective_name_coeff.end(); ++it) {
			int i = variable_names[it->first];
	        vector_c[i] = sense =='m' ? it->second : -it->second;
	    }

		// Solve Phase II LP.
		printf("*** PHASE II ***\n");
		printf("Problem data:\n");
		//print_tableau();
		simplex(basis);
	}
	
	print_elapsed_time();
		
	delete ipiv;
	ipiv = NULL;
	
	printf("Freeing memory.\n");
	delete [] matrix_A;
	delete [] vector_b;
	delete [] vector_c;
}

// This functions solves a system of linear equations of the form:
// A * x = b
// The result is stored in the _vector_bx, thus enough space must be allocated.
void LinearProgram::lineq_solve(double * _matrix_A, double * _vector_bx, double * _vector_cy, bool factorize = true)
{
	int info;
	int one_int = 1;
	if (ipiv == NULL) ipiv = new int[M];
	
	if (factorize) {
	    // LU factorization: A = P * L * U
	    // Arguments: num. of rows, num. of cols., matrix (on exit: LU factors), leading dimension of array, pivot indices, info.
	    dgetrf_(&M, &M, _matrix_A, &M, ipiv, &info);
	    if (info != 0) {
			// info < 0: argument "-info" has illegal value
	        // info > 0: diagonal element at "info" of the factor U from the factorization is exactly 0
	        printf("dgetrf_ failed with error code %d\n", (int) info);
	        exit(1);
	    }
			
		// Compute: A^T * y = c
	    // Arguments: transpose, order of matrix, num. of r.h.sides, LU factors in single matrix (from dgetrf),
	    // leading dimension of array, pivot indices from DGETRF, rhs (on exit: the solution), leading dimension of rhs, info.
		char transpose = 'T';
	    dgetrs_(&transpose, &M, &one_int, _matrix_A, &M, ipiv, _vector_cy, &M, &info);
	    if (info != 0) {
			// info < 0: argument "-info" has illegal value
	        printf("dgetrs_ #1 failed with error code %d\n", (int) info);
	        exit(1);
	    }
	}
	
    // Compute: A * x = b
	char transpose = 'N';
    dgetrs_(&transpose, &M, &one_int, _matrix_A, &M, ipiv, _vector_bx, &M, &info);
    if (info != 0) {
		// info < 0: argument "-info" has illegal value
        printf("dgetrs_ #1 failed with error code %d\n", (int) info);
        exit(1);
    }
}

// The (Revised) Simplex Algorithm.
// arg_basis : on input: initial basis; on result: final basis
int LinearProgram::simplex(set<int> & arg_basis)
{
	double * matrix_A_B = new double[M * M];
	double * matrix_A_N = new double[M * (N - M)];
	double * vector_c_B = new double[M];
	double * vector_c_N = new double[N - M];
	double * vector_bx = new double[M];
	double * vector_cy = new double[N];
	
	vector<int> basis, non_basis;
	for (int i = 0; i < N; i++) {
		if (arg_basis.count(i) == 0) non_basis.push_back(i); else basis.push_back(i);
	}
	
	unsigned long iteration_count = 0L;
	unsigned long iteration_limit = 1000000000L;
	while (iteration_count < iteration_limit) {
		iteration_count++;
		printf("\nSIMPLEX iteration: %d\n", iteration_count);
		
		// Print basis indices.
		// printf("Basis (%d): %d", basis.size(), basis[0]);
		// for (int i = 1; i < basis.size(); i++) {
		// 	printf(", %d", basis[i]);
		// }
		// printf("\n\n");
		
		// Split the matrix A into basis matrix (A_B) and non-basis matrix (A_N).
		int jB = 0, jN = 0;
		for (vector<int>::iterator it = basis.begin(); it != basis.end(); ++it) {
			memcpy(matrix_A_B + jB * M, matrix_A + (*it) * M, M * sizeof(double));
			vector_c_B[jB] = vector_c[*it];
			jB++;
		}
		for (vector<int>::iterator it = non_basis.begin(); it != non_basis.end(); ++it) {
			memcpy(matrix_A_N + jN * M, matrix_A + (*it) * M, M * sizeof(double));
			vector_c_N[jN] = vector_c[*it];
			jN++;
		}
		
		// printf("A_B=\n");
		// for (int i = 0; i < M; i++) {
		//     for (int j = 0; j < M; j++) {
		//         printf("%f\t", matrix_A_B[j * M + i]);
		//     }
		//     printf("\n");
		// }
		// printf("A_N=\n");
		// for (int i = 0; i < M; i++) {
		//     for (int j = 0; j < N-M; j++) {
		//         printf("%f\t", matrix_A_N[j * M + i]);
		//     }
		//     printf("\n");
		// }
		// printf("c_B=\n");
		// for (int i = 0; i < M; i++) {
		//    printf("%f\t", vector_c_B[i]);
		// }
		// printf("\n");
		// printf("c_N=\n");
		// for (int i = 0; i < N - M; i++) {
		//    printf("%f\t", vector_c_N[i]);
		// }
		// printf("\n");
	
		// Compute x = A_B^{-1} * b and y = (A_B^T)^{-1} * c_B
		memcpy(vector_bx, vector_b, M * sizeof(double));
		memcpy(vector_cy, vector_c_B, M * sizeof(double));
		lineq_solve(matrix_A_B, vector_bx, vector_cy);
	
		printf("solved x=\n");
		for (int i = 0; i < M; i++) {
		   printf("%f\t", vector_bx[i]);
		}
		printf("\nsolved y=\n");
		for (int i = 0; i < M; i++) {
		   printf("%f\t", vector_cy[i]);
			    }
			    printf("\n");
	
	    // Pricing: s = cN - (A_N)^T * y        
	    // Multiply matrix by vector: y = alpha * A * x + beta * y.
	    // Arguments: storage order, transpose, num. of rows, num. of cols., alpha, matrix A, l.d.a. of A, vect. x, incx, beta, y, incy.
		// NOTE 1: the first argument "storage order" is not in the original CBLAS specification
		// NOTE 2: the second argument should be 'T' according to the original CBLAS specification
	    cblas_dgemv(CblasColMajor, CblasTrans, M, N - M, -1.0, matrix_A_N, M, vector_cy, 1, 1.0, vector_c_N, 1);
	
		// Now vector_c_N contains the result s.
	
		printf("solved s=\n");
		for (int i = 0; i < N - M; i++) {
		   printf("%f\t", vector_c_N[i]);
		}
		printf("\n");
	
		// Bland's rule for pivoting.
		// 1) Choose entering index to be lexicographically first with s_j < 0.
		int entering_index = -1;
		for (int i = 0; i < N - M; i++) {
			if (vector_c_N[i] < 0.0) {
				entering_index = i;
				// Copy the entering column A_q into vector_cy.
				memcpy(vector_cy, matrix_A + (non_basis[i] * M), M * sizeof(double));
				//printf("Entering index=%d\n", non_basis[i]);
				break;
			}
		}

		if (entering_index == -1) {
			printf("OPTIMAL X FOUND.\n");
			double cost = 0.0;
            solution.clear();
            
			for (int i = 0; i < basis.size(); i++) {	
											
				// Applying shifts.		
				for (map<string, double>::iterator kt = var_shifts.begin(); kt != var_shifts.end(); ++kt) {
					if (variable_names[kt->first] == basis[i]) {
						printf("SHIFTING: %i by %f\n", i, kt->second);
						vector_bx[i] += kt->second;
					}
			
				}
				printf("var %i = %f\n", basis[i], vector_bx[i]);
                solution[variable_names_rev[basis[i]]] = vector_bx[i];
				cost += vector_c_B[i] * vector_bx[i];
			}
			for (int i = 0; i < non_basis.size(); i++) {
				double q = 0.0;
				// Applying shifts.		
				for (map<string, double>::iterator kt = var_shifts.begin(); kt != var_shifts.end(); ++kt) {
					if (variable_names[kt->first] == non_basis[i]) {
						printf("SHIFTING: var %i by %f\n", non_basis[i], kt->second);
						q += kt->second;
					}
				}
				printf("var %i = %f\n", non_basis[i], q);
                solution[variable_names_rev[non_basis[i]]] = q;
				cost += vector_c[non_basis[i]] * q;
			}
			objective_value = (sense == 'm') ? cost : -cost;
			printf("Value = %f\n", objective_value);
			
			arg_basis.clear();
			for (int i = 0; i < M; i++) arg_basis.insert(basis[i]);
			break;
		}
		
		// 2) Choose leaving index to be lexicographically first with min{ x_i / d_i, d_i > 0 }, 
		// d = A_B^{-1} * A(entering_index)
		
		// Solve: A_B * d = A(entering_index). Note that matrix_A_B already contains LU factors.
		lineq_solve(matrix_A_B, vector_cy, NULL, false);

		printf("solved d=\n");
		for (int i = 0; i < M; i++) {
		   printf("%f\t", vector_cy[i]);
		}
		printf("\n");
		
		int min_i = -1;
		double min_lbd = numeric_limits<double>::max();
		for (int i = 0; i < M; i++) {
			if (vector_cy[i] > 0.0) {
				double lbd = vector_bx[i] / vector_cy[i];
				if (lbd < min_lbd) {
					min_i = i;
					min_lbd = lbd;
				}
			}
		}
		
		if (min_i == -1) {
			printf("LP IS UNBOUNDED.\n");
			break;
		}
		
		int leaving_index = basis[min_i];
		//printf("Leaving index=%d\n", leaving_index);
		
		// Update basis.
		int tmp = basis[min_i];
		basis[min_i] = non_basis[entering_index];
		non_basis[entering_index] = tmp;
	}
	if (iteration_count == iteration_limit) {
		printf("MAXIMUM NUMER OF ITERATIONS REACHED.\n");
	}
	
	delete [] matrix_A_B;
	delete [] matrix_A_N;
	delete [] vector_c_B;
	delete [] vector_c_N;
	delete [] vector_bx;
	delete [] vector_cy;
	
    return 0;
}

// Write the solution to text file.
void LinearProgram::write(const string & filename)
{
    printf("Writing solution...\n");
    
    ofstream outfile;
    outfile.open(filename);
    
    outfile << "<?xml version = \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << endl;
    outfile << "<TestSolution><header objectiveValue=\"" << objective_value << "\"/>" << endl;
    
    for (map<string, double>::iterator it = solution.begin(); it != solution.end(); it++) {
        outfile << "<variable name=\"" << it->first.c_str() << "\" value=\"" << it->second << "\"/>" << endl;
    }
    outfile << "</TestSolution>" << endl;
    outfile.close();
}
