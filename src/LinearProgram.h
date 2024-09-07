// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#ifndef LINEAR_PROGRAM_H
#define LINEAR_PROGRAM_H

#include "Basis.h"

#include "Eigen/Dense"

#include <string_view>
#include <map>
#include <optional>
#include <vector>

namespace simplex
{
    struct Constraint
    {
        char type;
        double rhs;
        std::map<std::string, double> name_to_coeff;  // variable names and coefficients a_{ij}

        Constraint(char _type = '<', double _rhs = 0.0) : type(_type), rhs(_rhs) {}

        Constraint & add_term(std::string_view name, double a = 1.0);
        bool has_variable(std::string_view name) const;
        std::optional<double> get_coefficient(std::string_view name) const;
        Constraint & remove_term(std::string_view name);
        void negate_sides();
    };

    //
    // A standard form LP:
    // min c^T x
    // s.t.: Ax = b, 0 <= x <= UB
    //
    struct LinearProgram
    {
        friend class Presolve;
        friend class InitialBasis;
        friend class Simplex;

        void add_constraint(const std::string & name, Constraint && constraint);
        [[nodiscard]] const Constraint & get_constraint(const std::string & name) const;
        Eigen::Index add_variable(const std::string & var_name, double coeff = 0.0);
        void remove_variable(Eigen::Index var_id);
        bool has_variable(const std::string & var_name) const;
        Eigen::Index get_variable_id(const std::string & var_name) const;
        void set_sense(const char s) { sense = s; }
        void set_objective_label(const std::string & label) { objective_label = label; }
        void set_lower_bound(Eigen::Index var_id, double low_value);
        void set_upper_bound(Eigen::Index var_id, double low_value);

        /// @brief Copies the data read by parser into internal data structures.
        void initialize_tableau();

        /// @brief Print internal data structures on screen (for debugging).
        void print_tableau() const;

        /// @brief Add slack variables to inequality constraints, replacing them with equality.
        ///        The added slack variables are inserted into the basis.
        /// @return true if all the problem's constraints are inequalities
        void add_slack_variables_for_inequality_constraints();

        /// @brief Add artificial variables for each equality constraint, and return basis consisting of their indices.
        int add_artificial_variables_for_first_phase(Basis & basis);

        size_t get_num_vars() const { return variable_name_to_id.size(); }
        size_t get_num_rows() const { return constraints.size(); }

        /// @brief Returns the name of i-th artificial variable.
        static std::string get_artificial_variable(int i);

        /// @brief Apply upper-bound substitution to a given variable.
        void upper_bound_substitution(Eigen::Index var_id, double ub);
        bool has_non_trivial_upper_bounds() const { return m_has_UBS; }

        /// @brief Performs scaling of the matrix A and the vectors b and c.
        ///        Note that the original problem data in constraints and objective_coeff is not affected.
        void rescale_matrix();

    private:

        std::map<std::string, Constraint> constraints;
        std::map<Eigen::Index, double> objective_coeff;  // maps variable indices x_j to their c_j
        std::map<Eigen::Index, double> var_lbnd;         // LBs are set to 0 by variable substitutions in presovle
        std::map<Eigen::Index, double> var_ubnd;         // UBs are modified by presolve

        std::string objective_label;

        char sense = 'm';  // M=maximize / m=minimize

        Eigen::MatrixXd matrix_A;
        Eigen::VectorXd vector_b;
        Eigen::VectorXd vector_c;

        Eigen::Index next_variable_id = 0;
        std::map<std::string, Eigen::Index> variable_name_to_id;
        std::map<Eigen::Index, std::string> variable_id_to_name;

        std::vector<bool> ub_substitutions;
        std::vector<double> m_column_scaling_factors;
        std::vector<double> m_row_scaling_factors;
        std::map<std::string, Eigen::Index> m_slacks;
        std::map<std::string, Eigen::Index> m_artificials;

        double obj_value_shift = 0.0;  // constant term in objective function
        std::map<Eigen::Index, double> var_shifts;
        bool m_has_UBS = false;
    };
}  // namespace simplex

#endif
