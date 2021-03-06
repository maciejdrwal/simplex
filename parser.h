// Copyright (C) 2017-2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#ifndef _PARSER_H_
#define _PARSER_H_

class LinearProgram;

enum { LP_FILE, MPS_FILE } file_format;

bool run_parser_lp(const std::string& input, LinearProgram& lp);

// Old hand-written state machine parsers.
int parse_input_line_lp(const char * buffer, LinearProgram * lp);
int parse_input_line_mps(const char * buffer, LinearProgram * lp);

#endif
