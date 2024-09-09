// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#ifndef INITIAL_BASIS_H
#define INITIAL_BASIS_H

#include "LinearProgram.h"
#include "Basis.h"
#include "Simplex.h"

namespace simplex
{
    class InitialBasis
    {
    public:
        InitialBasis(LinearProgram & lp, Simplex & simplex) : m_lp(lp), m_simplex(simplex) {}

        Basis get_basis();

    private:

        LinearProgram & m_lp;
        Simplex & m_simplex;
        std::map<std::string, Eigen::Index> m_artificials;

        /// @brief Add artificial variables for each equality constraint and create an initial basis.
        ///        The problem must have only equality constraints at this point.
        ///        Available slack variables will be used instead of artificial variables when possible.
        /// @param basis The initial basis consisting of slack and/or artificial variables.
        /// @return Number of artificial variables added.
        size_t add_artificial_variables_for_initial_basis(Basis & basis);

    };
}  // namespace simplex

#endif
