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
        InitialBasis(LinearProgram & lp) : m_lp(lp) {}

        Basis get_basis();

    private:

        LinearProgram & m_lp;
    };
}  // namespace simplex

#endif
