// Copyright (C) 2017-2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#ifndef _PARSER_H_
#define _PARSER_H_

namespace simplex
{
    class LinearProgram;

    enum class FileFormat
    {
        LP_FILE,
        MPS_FILE
    };

    bool run_parser_lp(const std::string & input, LinearProgram & lp);
}

#endif
