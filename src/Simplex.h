// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "LinearProgram.h"

#include <map>

namespace simplex
{
    class Simplex
    {
    public:
        Simplex(LinearProgram & lp) : m_lp(lp) {}

        void solve();
        void write(const std::string & filename) const;

    private:

        LinearProgram & m_lp;

        std::map<std::string, double> m_solution;

        double m_objective_value = 0.0;

        int simplex(Basis & arg_basis);

        int select_entering_variable(const Eigen::VectorXd & s) const;
        int select_entering_variable_Bland(const Eigen::VectorXd & s) const;
        int select_entering_variable_most_neg(const Eigen::VectorXd & s) const;

        int select_leaving_variable(const Eigen::VectorXd & x, const Eigen::VectorXd & d,
                                    const Basis & basis, const Basis & non_basis,
                                    int entering_index) const;
        int select_leaving_variable_Bland(const Eigen::VectorXd & x, const Eigen::VectorXd & d) const;
        int select_leaving_variable_SUB(const Eigen::VectorXd & x, const Eigen::VectorXd & d, const Basis & basis,
                                        const Basis & non_basis, int entering_index) const;

        void solution_found(Eigen::VectorXd vector_bx, const Eigen::VectorXd & vector_c_B, const Basis & basis,
                            const Basis & non_basis);
    };

}  // namespace simplex

#endif
