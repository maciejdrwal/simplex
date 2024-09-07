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

namespace simplex
{
    class InitialBasis
    {
    public:
        InitialBasis(LinearProgram & lp, Simplex & simplex) : m_lp(lp), m_simplex(simplex) {}

        Basis get_basis() const;

    private:

        LinearProgram & m_lp;
        Simplex & m_simplex;
    };
}  // namespace simplex

#endif
