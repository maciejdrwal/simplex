// Copyright (C) 2024 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#ifndef PARSER_H
#define PARSER_H

#include "Simplex.h"

namespace simplex
{
    bool run_parser_lp(std::string_view input, simplex::LinearProgram & lp);
}

#endif
