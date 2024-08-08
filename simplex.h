// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <map>
#include <optional>
#include <set>
#include <string_view>
#include <vector>

#include "Eigen/Dense"

namespace simplex
{
    struct Constraint
    {
        char type;
        double rhs;
        std::map<std::string, double> name_to_coeff;  // variable names and coefficients a_{ij}

        Constraint(char _type = '<', double _rhs = 0.0) : type(_type), rhs(_rhs) {}

        void add_term(std::string_view name, double a)
        {
            const auto it = name_to_coeff.find(name.data());
            if (it == name_to_coeff.end())
            {
                name_to_coeff[name.data()] = a;  // add new term
            }
            else
            {
                it->second += a;  // add to existing term
            }
        }

        bool has_variable(std::string_view name) const { return name_to_coeff.count(name.data()) > 0; }

        std::optional<double> get_coefficient(std::string_view name) const
        {
            const auto it = name_to_coeff.find(name.data());
            return it == name_to_coeff.end() ? std::nullopt : std::optional(it->second);
        }

        void remove_term(std::string_view name)
        {
            const auto it = name_to_coeff.find(name.data());
            if (it != name_to_coeff.end())
            {
                name_to_coeff.erase(it);
            }
        }
    };

    //
    // A standard form LP:
    // min c^T x
    // s.t.: Ax = b, 0 <= x <= UB
    //
    struct LinearProgram
    {
        friend class Presolve;

        std::map<std::string, Constraint> constraints;
        std::map<Eigen::Index, double> objective_coeff;  // maps variable indices x_j to their c_j
        std::map<Eigen::Index, double> var_lbnd;         // LBs are set to 0 by variable substitutions in presovle
        std::map<Eigen::Index, double> var_ubnd;         // UBs are modified by presolve
        double obj_value_shift;                          // constant term in objective function

        using Basis = std::vector<Eigen::Index>;

        LinearProgram() : obj_value_shift(0.0), sense('m'), objective_value(0.0) {}

        /// @brief Add slack variables to inequality constraints, replacing them with equality.
        ///        The added slack variables are inserted into the basis.
        /// @return true if all the problem's constraints are inequalities
        bool add_slack_variables_for_inequality_constraints(Basis & basis);

        /// @brief Add artificial variables for each equality constraint, and return basis consisting of their indices.
        int add_artificial_variables_for_first_phase(Basis & basis);

        void solve();
        void write(const std::string & filename) const;

        Eigen::Index add_variable(const std::string & var_name, double coeff = 0.0);
        void remove_variable(Eigen::Index var_id);
        bool has_variable(const std::string & var_name) const;
        void set_sense(const char s) { sense = s; }
        void set_objective_label(const std::string & label) { objective_label = label; }

        std::string objective_label;

    private:
        Eigen::MatrixXd matrix_A;
        Eigen::VectorXd vector_b;
        Eigen::VectorXd vector_c;

        Eigen::Index next_variable_id = 0;

        char sense;  // M=maximize / m=minimize
        double objective_value;

        std::vector<bool> ub_substitutions;

        std::map<std::string, Eigen::Index> variable_name_to_id;
        std::map<Eigen::Index, std::string> variable_id_to_name;
        std::map<Eigen::Index, double> var_shifts;
        std::map<std::string, double> solution;

        void initialize_tableau();
        void print_tableau() const;
        int simplex(Basis & arg_basis);
        void lineq_solve(Eigen::MatrixXd & matrix_A, Eigen::VectorXd & vector_bx, Eigen::VectorXd & vector_cy,
                         bool factorize);
        int select_entering_variable_Bland(const Eigen::VectorXd & vector_c_N) const;
        int select_entering_variable_most_neg(const Eigen::VectorXd & vector_c_N) const;
        int select_leaving_variable_Bland(Eigen::VectorXd & vector_bx, Eigen::VectorXd & vector_cy,
                                          Eigen::MatrixXd & matrix_A_B);
        int select_leaving_variable_SUB(double * vector_bx, double * vector_cy, double * matrix_A_B,
                                        int entering_index);
        void upper_bound_substitution(Eigen::Index var_id, double ub);
        void solution_found(Eigen::VectorXd & vector_bx, const Eigen::VectorXd & vector_c_B, 
                            const Basis & basis, const Basis & non_basis);
    };
}  // namespace simplex

#endif
