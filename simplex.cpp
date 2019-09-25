// Copyright (C) 2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#include <chrono>
#include <csignal>
#include <fstream>
#include <iostream>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#endif

#ifdef __linux__
#include <clbas.h>
#include <f77blas.h>
#endif

#include "presolve.h"
#include "simplex.h"
#include "utils.h"

using namespace std;

int * LinearProgram::ipiv = NULL;

chrono::time_point<chrono::steady_clock> t_start, t_end;

// Print internal data structures on screen (for debugging).
void LinearProgram::print_tableau() const
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

void print_matrix(double * matrix, int n_rows, int m_cols)
{
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            cout << scientific << matrix[j * n_rows + i] << '\t';
        }
        cout << fixed << endl;
    }
}

template <typename T>
void print_vector(T * vect, int n)
{
    for (int i = 0; i < n; i++) {
        cout << scientific << vect[i] << '\t';
    }
    cout << fixed << endl;
}

// Copies the data read by parser into internal data structures: A, b.
void LinearProgram::initialize_tableau()
{
    int j = 0;
    for (map<string, shared_ptr<Constraint> >::iterator it = constraints.begin(); it != constraints.end(); ++it) {
        for (map<string, double>::iterator jt = it->second->name_coeff.begin(); jt != it->second->name_coeff.end(); ++jt) {
            int i = variable_names_to_ids[jt->first] * M + j;
            matrix_A[i] = jt->second;
        }
        vector_b[j++] = it->second->rhs;
    }
    objective_value = 0.0;
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
// Initializes internal data structures and runs 2 phases of simplex algorithm.
void LinearProgram::solve()
{        
    measure_time_start();
    
    presolve(this);
    
    M = constraints.size();
    N = variable_names_to_ids.size();
    
    printf("problem size:\nN=%d, M=%d\n", N, M);
    
    cout << "Allocating memory..." << endl;
    
    vector_b = new double[M];
    vector_c = new double[N + 2 * M];   // add reserve for slack/artificial variables
    matrix_A = new double[M * (N + 2 * M)];
    
    memset(vector_b, 0, M * sizeof(double));
    memset(vector_c, 0, (N + 2 * M) * sizeof(double));
    memset(matrix_A, 0, (M * (N + 2 * M)) * sizeof(double));
                
    set<int> init_basis;     // Initial basis.
        
    // Add slack variables to inequality constraints, replacing them with equality.
    int num_of_slack_vars = 0;
    for (map<string, shared_ptr<Constraint> >::iterator it = constraints.begin(); it != constraints.end(); ++it) {
        if (it->second->type != '=') {
            string var_name("__SLACK" + tostr<int>(num_of_slack_vars++));
            if (it->second->type == '<') {
                it->second->name_coeff[var_name] = 1.0;
            }
            if (it->second->type == '>') {
                it->second->name_coeff[var_name] = -1.0;
            }
            it->second->type = '=';
            add_variable(var_name);
            
            // If all constraints were inequalities then use slacks for initial basis.
            if (all_inequalities) init_basis.insert(variable_names_to_ids[var_name]);
            printf("added slack variable:%s\n", var_name.c_str());
        }
    }
    int num_of_primal_vars = objective_name_coeff.size();
    
    N = variable_names_to_ids.size();   // this might have been modified by adding slack variables
    
    if (!all_inequalities) {
        // Construct artificial variables for Phase I.
        int i = 0;
        printf("Not all constrains are inequalities A <= b, adding artificial variables.\n");
        for (map<string, shared_ptr<Constraint> >::iterator it = constraints.begin(); it != constraints.end(); ++it) {
            if (it->second->rhs < 0.0) {
                it->second->rhs = -it->second->rhs;
                for (map<string, double>::iterator jt = it->second->name_coeff.begin(); jt != it->second->name_coeff.end(); ++jt) {
                    jt->second = -jt->second;
                }
                // Note: at this point there should be no inequality constraints
                if (it->second->type == '<' || it->second->type == '>') {
                    throw "All constraints must be equality at this point.";
                }
            }
            string var_name("__ARTIFICIAL" + tostr<int>(i++));
            it->second->name_coeff[var_name] = 1.0;
            add_variable(var_name);
            init_basis.insert(variable_names_to_ids[var_name]);
            printf("added artificial variable: %s\n", var_name.c_str());
        }
        // Objective function: sum of artificial variables.
        for (int i = 0; i < M; i++) vector_c[num_of_primal_vars + num_of_slack_vars + i] = 1.0;
        N = variable_names_to_ids.size();   // this has changed after adding artificial variables
        
        // Solve Phase I LP.
        printf("*** PHASE I ***\n");
        initialize_tableau();
        print_tableau();   
        simplex(init_basis);
        
        printf("obj=%f\n", objective_value);
        
        // Remove artificial variables.
        i = 0;
        for (map<string, shared_ptr<Constraint> >::iterator it = constraints.begin(); it != constraints.end(); ++it) {
            string var_name = string("__ARTIFICIAL") + tostr<int>(i++);
            int var_id = variable_names_to_ids[var_name];

            // Check if some artificial variable remained in the basis.
            if (init_basis.count(var_id) > 0) {
                // TODO: handle this case
                
                init_basis.erase(variable_names_to_ids[var_name]);
                for (map<string, int>::iterator jt = variable_names_to_ids.begin(); jt != variable_names_to_ids.end(); ++jt) {
                    if (init_basis.count(jt->second) == 0 && (strncmp(jt->first.c_str(), "__ARTIFICIAL", 12) != 0)) {
                        init_basis.insert(jt->second);
                        printf("replaced variable %d (%s) by %d (%s) in the basis\n",  var_id, var_name.c_str(), jt->second, jt->first.c_str());
                        break;
                    }
                }
                
                printf("Artificial variable %s in the basis. Aborting.\n", var_name.c_str());
                exit(1);
            }
            
            it->second->name_coeff.erase(var_name);
            variable_names_to_ids.erase(var_name);
            variable_ids_to_names.erase(var_id);
            printf("removed artificial variable:%s\n", var_name.c_str());
        }
        N = variable_names_to_ids.size();   // this has changed after removing artificial variables
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
            int i = variable_names_to_ids[it->first];
            vector_c[i] = (sense =='m') ? it->second : -it->second;
        }

        // Solve Phase II LP.
        printf("*** PHASE II ***\n");
        //printf("Problem data:\n");
        print_tableau();
        simplex(init_basis);
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
// The result is stored in _vector_bx, thus enough space must be allocated.
// If factorize = true then also the following system is solved:
// A^T * y = c
// The result is stored in _vector_cy, thus enough space must be allocated.
// If factorize = false then _matrix_A must contain LU factors of A on input.
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
        printf("dgetrs_ #2 failed with error code %d\n", (int) info);
        exit(1);
    }
}

// Bland's rule for pivoting.
// Choose entering index to be lexicographically first with s_j < 0.
int LinearProgram::select_entering_variable_Bland(double * vector_c_N)
{
    for (int i = 0; i < N - M; i++) {
        if (vector_c_N[i] < 0.0) return i;
    }
    return -1;
}

// Bland's rule for pivoting.
// Choose leaving index to be lexicographically first with min{ x_i / d_i, d_i > 0 }, 
// d = A_B^{-1} * A(entering_index)
int LinearProgram::select_leaving_variable_Bland(double * vector_bx, double * vector_cy, double * matrix_A_B)
{
    // Solve: A_B * d = A(entering_index). Note that matrix_A_B already contains LU factors.
    lineq_solve(matrix_A_B, vector_cy, NULL, false);

    // cout << "solved d=\n";
    // print_vector<double>(vector_cy, M);
    
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
    return min_i;
}

// Choose the most negative value among s_j < 0.
int LinearProgram::select_entering_variable_most_neg(double * vector_c_N)
{
    int min_i = -1;
    double min_value = numeric_limits<double>::max();
    for (int i = 0; i < N - M; i++) {
        if (vector_c_N[i] < 0.0 && vector_c_N[i] < min_value) {
            min_i = i;
            min_value = vector_c_N[i];
        } 
    }
    return min_i;
}

int LinearProgram::select_leaving_variable_GUB(double * vector_bx, double * vector_cy, double * matrix_A_B, int entering_index)
{
    // Solve: A_B * d = A(entering_index). Note that matrix_A_B already contains LU factors.
    // Result d is stored in vector_cy.
    lineq_solve(matrix_A_B, vector_cy, NULL, false);

    bool need_substitution = false;
    int min_i = -1;
    double min_theta = numeric_limits<double>::max();
    for (int i = 0; i < M; i++) {
        if (vector_cy[i] > 0.0) {
            double t = vector_bx[i] / vector_cy[i];
            if (t < min_theta) {
                min_i     = i;
                min_theta = t;
            }
        }
    }

    for (int i = 0; i < M; i++) {
        if (vector_cy[i] < 0.0) {
            double ub = var_ubnd[variable_ids_to_names[basis[i]]];
            double t  = (ub - vector_bx[i]) / (-vector_cy[i]);
            if (t < min_theta) {
                min_i     = i;
                min_theta = t;
                need_substitution = true;
            }
        }
    }

    double ub_s = var_ubnd[variable_ids_to_names[non_basis[entering_index]]];

    if (ub_s < min_theta) {
        upper_bound_substitution(non_basis[entering_index], ub_s);
        return -2;  // do not perform pivot
    }
    else {
        if (need_substitution) {
            int k = basis[min_i];
            upper_bound_substitution(k, var_ubnd[variable_ids_to_names[k]]);
        }
    }

    return min_i;   // perform usual pivot (except if min_i == -1)
}

void LinearProgram::upper_bound_substitution(int var_id, double ub)
{
    // Modify the objective function: c_s x_s -> -c_s x_s
    vector_c[var_id] = -vector_c[var_id];

    // Modify constraints: a_{js} x_s -> -a_{js} x_s, b_j -> b_j - a_{js}*u_s
    for (int i = 0; i < M; i++) {
        double * a_ptr = matrix_A + var_id * M + i;
        double   a     = *a_ptr;
        vector_b[i] = vector_b[i] - a * ub;
        *a_ptr = -a;
    }
    
    ub_substitutions[var_id] = !ub_substitutions[var_id];
}

void LinearProgram::solution_found(double * vector_bx, double * vector_cy, double * vector_c_B)
{
    double cost = 0.0;
    solution.clear();

    for (int i = 0; i < basis.size(); i++) {    
        int var_id = basis[i];
        string var_name = variable_ids_to_names[var_id];

        // Reverse UB-substitution
        if (ub_substitutions[var_id]) {
            cout << "REVERSING UB-subst of " << var_name << endl;
            upper_bound_substitution(var_id, var_ubnd[var_name]);
        }

        // Applying shifts
        for (map<string, double>::iterator kt = var_shifts.begin(); kt != var_shifts.end(); ++kt) {
            if (variable_names_to_ids[kt->first] == var_id) {
                cout << "SHIFTING: basis var " << var_name << " by " << kt->second << endl;
                vector_bx[i] += kt->second;
            }
        }
        
        cout << "var (basis) " << var_name << " = " << vector_bx[i] << endl;
        solution[var_name] = vector_bx[i];
        cost += vector_c_B[i] * vector_bx[i];
    }
    for (int i = 0; i < non_basis.size(); i++) {
        double q = 0.0;
        int var_id = non_basis[i];
        string var_name = variable_ids_to_names[var_id];

        // Reverse UB-substitution
        if (ub_substitutions[var_id]) {
            double ub = var_ubnd[var_name];
            upper_bound_substitution(var_id, ub);
            q += ub;
            cout << "REVERSING UB-subst of " << var_name << " ub= " << ub << endl;
        }

        // Applying shifts.     
        for (map<string, double>::iterator kt = var_shifts.begin(); kt != var_shifts.end(); ++kt) {
            if (variable_names_to_ids[kt->first] == var_id) {
                cout << "SHIFTING: non-basis var " << var_name << " by " <<  kt->second << endl;
                q += kt->second;
            }
        }
        
        cout << "var (non-basis) " << var_name << " = " << q << endl;
        solution[var_name] = q;
        cost += vector_c[non_basis[i]] * q;
    }
    objective_value = (sense == 'm') ? cost : -cost;
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
    
    basis.clear();
    non_basis.clear();
    for (int i = 0; i < N; i++) {
        if (arg_basis.count(i) == 0) non_basis.push_back(i); else basis.push_back(i);
    }

    // Initialize upper-bound substitutions.
    ub_substitutions = vector<bool>(N, false);
    
    unsigned long iteration_count = 0L;
    unsigned long iteration_limit = 1000000000L;
    while (iteration_count < iteration_limit) {
        iteration_count++;
        cout << "\nSIMPLEX iteration: " << iteration_count << endl;
        cout << "Basis: ";
        print_vector<int>(basis.data(), basis.size());
        
        // Split the matrix A into basis matrix (A_B) and non-basis matrix (A_N).
        // TODO: no need to copy whole vectors, as only 1 element has changed in basis.
        // However note that matrix_A, vector_b and vector_c might have changed
        // due to applying upper-bound substitutions.
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
        
        cout << "A_B=\n";
        print_matrix(matrix_A_B, M, M);
        cout << "A_N=\n";
        print_matrix(matrix_A_N, M, N-M);
        cout << "c_B=\n";
        print_vector<double>(vector_c_B, M);
        cout << "c_N=\n";
        print_vector<double>(vector_c_N, N - M);
    
        // Compute x = A_B^{-1} * b and y = (A_B^T)^{-1} * c_B
        memcpy(vector_bx, vector_b, M * sizeof(double));
        memcpy(vector_cy, vector_c_B, M * sizeof(double));
        lineq_solve(matrix_A_B, vector_bx, vector_cy);
    
        cout << "solved x=\n";
        print_vector<double>(vector_bx, M);
        cout << "solved y=\n";
        print_vector<double>(vector_cy, M);
    
        // Pricing (reduced costs): s = c_N - (A_N)^T * y        
        // Multiply matrix by vector: y = alpha * A * x + beta * y.
        // Arguments: storage order, transpose, num. of rows, num. of cols., alpha, matrix A, l.d.a. of A, vect. x, incx, beta, y, incy.
        // NOTE 1: the first argument "storage order" is not in the original CBLAS specification
        // NOTE 2: the second argument should be 'T' according to the original CBLAS specification
        cblas_dgemv(CblasColMajor, CblasTrans, M, N - M, -1.0, matrix_A_N, M, vector_cy, 1, 1.0, vector_c_N, 1);
    
        // Now vector_c_N contains the result s.
    
        cout << "reduced costs s=\n";
        print_vector<double>(vector_c_N, N - M);
        
        // 1) Select the entering variable.
        int entering_index = select_entering_variable_Bland(vector_c_N);
        //int entering_index = select_entering_variable_most_neg(vector_c_N);
        
        if (entering_index == -1) {
            solution_found(vector_bx, vector_cy, vector_c_B);
            arg_basis.clear();
            for (int i = 0; i < M; i++) arg_basis.insert(basis[i]);

            cout << "OPTIMAL X FOUND.\nValue = " << objective_value << endl;
            break;
        }

        // Copy the entering column A_q into vector_cy.
        memcpy(vector_cy, matrix_A + (non_basis[entering_index] * M), M * sizeof(double));
        
        // 2) Select the leaving variable.
        // int leaving_index = select_leaving_variable_Bland(vector_bx, vector_cy, matrix_A_B);
        int leaving_index = select_leaving_variable_GUB(vector_bx, vector_cy, matrix_A_B, entering_index);

        if (leaving_index == -1) {
            cout << "LP IS UNBOUNDED.\n";
            break;
        }

        if (leaving_index != -2) {
            // Update basis.
            int tmp = basis[leaving_index];
            basis[leaving_index] = non_basis[entering_index];
            non_basis[entering_index] = tmp;
        }
    }

    if (iteration_count == iteration_limit) {
        throw "MAXIMUM NUMER OF SIMPLEX ITERATIONS REACHED";
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
void LinearProgram::write(const string & filename) const
{
    cout << "Writing solution...\n";
    
    ofstream outfile;
    outfile.open(filename);
    
    outfile << "<?xml version = \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << endl;
    outfile << "<TestSolution>\n<header objectiveValue=\"" << objective_value << "\"/>" << endl;
    
    for (map<string, double>::const_iterator it = solution.begin(); it != solution.end(); it++) {
        outfile << "<variable name=\"" << it->first.c_str() << "\" value=\"" << it->second << "\"/>" << endl;
    }
    outfile << "</TestSolution>" << endl;
    outfile.close();
}

// Note: each variable has ID that corresponds to that variable's position
// in the computation array double * vector_c, and is also used, e.g., to 
// indicate basic/non-basic variables. Thus it is currently not allowed
// to add variables after removing some other variables.
// TODO: change this in order to allow problem manipulation after solving.
void LinearProgram::add_variable(const string & var_name)
{
    static int next_id = 0;
    variable_names_to_ids[var_name] = next_id;
    variable_ids_to_names[next_id] = var_name;
    next_id++;
    var_lbnd[var_name] = 0.0;   // set default bounds
    var_ubnd[var_name] = std::numeric_limits<double>::max();
}

bool LinearProgram::has_variable(const string & var_name) const
{
    return (0 != variable_names_to_ids.count(var_name));
}
