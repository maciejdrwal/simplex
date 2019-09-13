// Copyright (C) 2017-2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sstream>

#include "parser.h"
#include "utils.h"

#define DEBUG_MODE 1
#define debug_printf(fmt, ...) do { if (DEBUG_MODE) printf(fmt, ##__VA_ARGS__); } while (0)
#define strncpy0(src, dst, num) { strncpy(src, dst, num); src[num] = '\0'; }
#define _abs(x) ((x) < 0.0 ? -(x) : (x))
#define _isfloatzero(x) (_abs(x) < 1e-9 ? 1 : 0)

int parser_state = 0;
int bounds_count = 0;
string constr_name;
map<string, double> constr_name_coeff;
	
// note: this function removes endline character at the end
char * trimwhitespace(char *str) 
{
    char *end;
    while(isspace(*str)) str++;   // trim leading space
    if (*str == 0) { return str; }    // all spaces?
    end = str + strlen(str) - 1;    // trim trailing space
    while(end > str && isspace(*end)) { end--; }
    *(end + 1) = 0;     // write new null terminator
    return str;
}

void apply_shifts(LinearProgram * lp) {
	// printf("all variables:\n");
	// for (map<string, int>::iterator it = lp->variable_names.begin(); it != lp->variable_names.end(); ++it) {
	// 	printf("%s -> %i\n", it->first.c_str(), it->second);
	// }
	//
	
	for (map<string, double>::iterator it = lp->objective_name_coeff.begin(); it != lp->objective_name_coeff.end(); ++it) {
		
		// Applying shifts.
		for (map<string, double>::iterator kt = lp->var_shifts.begin(); kt != lp->var_shifts.end(); ++kt) {
			if ((kt->first == it->first) && (_isfloatzero(kt->second))) {
				printf("APPLYING ZERO SHIFT: %s\n", it->first.c_str());
				it->second = -it->second;
			}
		}
		
		printf("%f * %s + ", it->second, it->first.c_str());
	}
	printf("\n\nConstraints:\n");
	for (map<string, Constraint*>::iterator it = lp->constraints.begin(); it != lp->constraints.end(); ++it) {
		char _type = it->second->type;
		double _rhs = it->second->rhs;
		if ((_type == '<' && _rhs < 0.0) ||
			(_type == '>' && _rhs > 0.0) ||
			(_type == '=')) lp->all_inequalities = false;

		printf("%s : ", it->first.c_str());
		for (map<string, double>::iterator jt = it->second->name_coeff.begin(); jt != it->second->name_coeff.end(); ++jt) {
			printf("%f %s + ", jt->second, jt->first.c_str());
			
			// Applying shifts.		
			for (map<string, double>::iterator kt = lp->var_shifts.begin(); kt != lp->var_shifts.end(); ++kt) {
				if (kt->first == jt->first) {
					if (_isfloatzero(kt->second)) {
						printf("APPLYING ZERO SHIFT: %s -> %f\n", kt->first.c_str(), kt->second);
						jt->second = -jt->second;
					}
					else {
						printf("APPLYING POS. SHIFT: %s -> %f\n", kt->first.c_str(), kt->second);
						it->second->rhs -= kt->second * jt->second;
					}
				}		
			}
		}
		printf("%c %f\n", it->second->type, it->second->rhs);
	}
}


// This function implements the parser of .lp files.
// The format for storing LP/ILP is used by CPLEX solver.
// TODO: handle syntax errors, etc.
int parse_input_line_lp(const char * buffer, LinearProgram * lp) 
{
    char _buf[1024], *line_buf;
	line_buf = _buf;
    
    debug_printf("parser buffer: %s", buffer);
    
    // copy the line data to internal buffer
    strcpy(line_buf, buffer);
    line_buf = trimwhitespace(line_buf);
    if (strcmp(buffer, "\n") == 0 || strlen(line_buf) == 0) {
        debug_printf("parser: empty line\n");
        return 0;
    }
    
    if (line_buf[0] == '\\' || line_buf[0] == '#' || line_buf[0] == '/') {
        debug_printf("parser: comment line\n");
        return 0;
    }
    
    if (strncmp(line_buf, "Maximize", 8) == 0) {
        lp->sense = 'M';
		parser_state = 1;
        return 0;
    }
	else if (strncmp(line_buf, "Minimize", 8) == 0) {
        lp->sense = 'm';
		parser_state = 1;
        return 0;
    }
    else if (strncmp(line_buf, "Subject To", 10) == 0) {
        parser_state = 2;
        return 0;
    }
    else if (strncmp(line_buf, "Bounds", 6) == 0) {
        parser_state = 3;
        return 0;
    }
    else if (strncmp(line_buf, "End", 3) == 0) {
		
		debug_printf("Objective label: %s\n", lp->objective_label.c_str());
		debug_printf("Number of variables: %d\n", lp->variable_names.size());
		
		apply_shifts(lp);
		
		// Reset the internal parser state, to make it ready to be used on a different file.
		parser_state = 0;
		bounds_count = 0;
		
        return 0;
    }
    
	char label[128], contents[1024], term[128], coeff[128], var_name[128], lb[128], ub[128];
	const char *contents_ptr, *contents_next_ptr;
	int span;
	double sign = 1.0;
	
    switch (parser_state) {
        case 1:
		// handle the objective function
		span = strcspn(line_buf, ":");
		if (span == strlen(line_buf)) {
			span = -1;
		}
		else {
			strncpy0(label, line_buf, span);
			if (lp->objective_label.length() == 0) lp->objective_label = string(label);
		}
		strcpy(contents, trimwhitespace(line_buf + span + 1));
			
		// tokenize the contents of the line, of the form: coeff1 var_name1 + coeff2 var_name2 + ...
		contents_ptr = contents;
		while (contents_ptr != NULL) {
			sign = 1.0;
			if (*contents_ptr == '-') { contents_ptr++; sign = -1.0; }
			if (*contents_ptr == '+') contents_ptr++;
			
			contents_next_ptr = contents_ptr;
			contents_ptr = strpbrk(contents_next_ptr, "+-");
			if (contents_ptr == NULL) {
				strcpy(term, contents_next_ptr);
			}
			else {
				strncpy0(term, contents_next_ptr, contents_ptr - contents_next_ptr);
			}
			strcpy(term, trimwhitespace(term));
			if (strchr(term, ' ') == NULL) {
				strcpy(coeff, "1");
				strcpy(var_name, term);
			}
			else {
				strncpy0(coeff, term, strchr(term, ' ') - term);
				strcpy(var_name, strchr(term, ' '));
			}
			strcpy(var_name, trimwhitespace(var_name));
			
			debug_printf("term: %s, coeff: %f, var_name: %s\n", term, sign * atof(coeff), var_name);
			
			lp->variable_names[var_name] = lp->variable_names.size();
			if (strlen(coeff) == 0) {
				lp->objective_name_coeff[var_name] = sign;
			}
			else {
				lp->objective_name_coeff[var_name] = sign * atof(coeff);
			}
		}		
        return 0;
        
        case 2:
		// handle a constraint
		span = strcspn(line_buf, ":");
		if (span == strlen(line_buf)) {
			span = -1;
		}
		else {
			strncpy0(label, line_buf, span);
			if (constr_name.length() == 0) constr_name = string(label);
		}
				
		strcpy(contents, trimwhitespace(line_buf + span + 1));
		
		// tokenize the contents of the line
		contents_ptr = contents;
		while (contents_ptr != NULL) {
			sign = 1.0;
			if (*contents_ptr == '-') { contents_ptr++; sign = -1.0; }
			if (*contents_ptr == '+') contents_ptr++;
			if (*contents_ptr == '=' || *contents_ptr == '<' || *contents_ptr == '>') {
				// add new constraint information to the lp object
				strcpy(term, contents_ptr + strspn(contents_ptr, "<>="));
				
				Constraint * constraint_ptr = new Constraint(*contents_ptr, atof(term), constr_name_coeff);
				lp->constraints[constr_name] = constraint_ptr;
				
				// update set of defined variables if necessary
				for (map<string, double>::iterator it = constr_name_coeff.begin(); it != constr_name_coeff.end(); ++it) {
			        if (lp->variable_names.count(it->first) == 0) {
					    debug_printf("parser: adding new variable via constraint:%s\n", it->first.c_str());
						lp->variable_names[it->first] = lp->variable_names.size();
					}
				}
				
				if ((*contents_ptr == '<' && constraint_ptr->rhs < 0.0) ||
					(*contents_ptr == '>' && constraint_ptr->rhs > 0.0) ||
					(*contents_ptr == '=')) lp->all_inequalities = false;
				
				constr_name = string();
				constr_name_coeff.clear();				
				debug_printf("parser: new constraint added\n");
				return 0;
			}
			
			contents_next_ptr = contents_ptr;
			contents_ptr = strpbrk(contents_next_ptr, "+-<>=");
									
			if (contents_ptr == NULL) {
				strcpy(term, contents_next_ptr);
			}
			else {
				strncpy0(term, contents_next_ptr, contents_ptr - contents_next_ptr);
			}
			strcpy(term, trimwhitespace(term));	
			
			if (strchr(term, ' ') == NULL) {
				strcpy(coeff, "1");
				strcpy(var_name, term);
			}
			else {
				strncpy0(coeff, term, strchr(term, ' ') - term);
				strcpy(var_name, strchr(term, ' '));
			}
			strcpy(var_name, trimwhitespace(var_name));
			
			debug_printf("term: %s, coeff: %f, var_name: %s\n", term, sign * atof(coeff), var_name);
			
			constr_name_coeff[var_name] = sign * atof(coeff);
		}
        
        return 0;
        
        case 3:
		// handle bounds
		contents_next_ptr = line_buf;
		contents_ptr = strchr(contents_next_ptr, '<');
		if (contents_ptr != NULL) {
			strncpy0(lb, contents_next_ptr, contents_ptr - contents_next_ptr);
			contents_next_ptr = contents_ptr + strspn(contents_ptr, "<= ");
			contents_ptr = strchr(contents_next_ptr, '<');
			if (contents_ptr != NULL) {
				strncpy0(var_name, contents_next_ptr, contents_ptr - contents_next_ptr);
				strcpy(var_name, trimwhitespace(var_name));
				strcpy(ub, contents_ptr + strspn(contents_ptr, "<= "));
				
				double lbf = atof(lb), ubf = atof(ub);
				if (lbf < 0.0) {
					// negative lower bound: substitute z = x + lbf, z >= 0
					lp->var_shifts[var_name] = lbf;
				}
				else if (lbf > 0.0) {
					// add ordinary constraint
					constr_name_coeff[var_name] = 1.0;
					Constraint * constraint = new Constraint('>', lbf, constr_name_coeff);
					lp->constraints[string("#LBOUND") + tostr<int>(bounds_count++)] = constraint;
					constr_name_coeff.clear();
				}
				
				char _type = '<';
				if (ubf < 0.0) {
					// negative upper bound: substitute z = -x
					lp->var_shifts[var_name] = 0.0;
					ubf = -ubf;
					_type = '>';
				}
				//add ordinary constraint
				constr_name_coeff[var_name] = 1.0;
				Constraint * constraint = new Constraint(_type, ubf, constr_name_coeff);				
				lp->constraints[string("#UBOUND") + tostr<int>(bounds_count++)] = constraint;
				constr_name_coeff.clear();
			}
		}
		
		debug_printf("lb: %s, var: %s, ub: %s\n", lb, var_name, ub);
		
        return 0;        
    }
    
    return 0;
}

// This function implements the parser of .mps (Mathematical Programming System) files.
// TODO: handle syntax errors, etc.
int parse_input_line_mps(const char *buffer, LinearProgram *lp)
{
    char _buf[1024], *line_buf;
	line_buf = _buf;

    debug_printf("parser: buffer: %s", buffer);

    // copy the line data to internal buffer
    strcpy(line_buf, buffer);
    line_buf = trimwhitespace(line_buf);
    if (strcmp(buffer, "\n") == 0 || strlen(line_buf) == 0) {
        debug_printf("parser: empty line\n");
        return 0;
    }

    if (line_buf[0] == '*') {
        debug_printf("parser: comment line\n");
        return 0;
    }

    if (strncmp(line_buf, "ROWS", 4) == 0) {
        if (parser_state != 1) { parser_state = 1; return 0; }
    }
    else if (strncmp(line_buf, "COLUMNS", 7) == 0) {
        if (parser_state != 2) { parser_state = 2; return 0; }
    }
    else if (strncmp(line_buf, "RHS", 3) == 0) {
        if (parser_state != 3) { parser_state = 3; return 0; }
    }
    else if (strncmp(line_buf, "BOUNDS", 6) == 0) {
        if (parser_state != 4) { parser_state = 4; return 0; }
    }
    else if (strncmp(line_buf, "ENDATA", 6) == 0) {
		
		debug_printf("Objective label: %s\n", lp->objective_label.c_str());
		debug_printf("Number of variables: %d\n", lp->variable_names.size());
		
		apply_shifts(lp);
		
		// Reset the internal parser state, to make it ready to be used on a different file.
		parser_state = 0;
		bounds_count = 0;
		
        return 0;
    }

    char rhs_name[16], col_name[16], row_name1[16], data1[16], row_name2[16], data2[16], row_type[16];
    double f1, f2;
    int cidx, ridx;
    string cn, rn;

    switch (parser_state) {
        case 1:
        sscanf(line_buf, "%16s %16s", &row_type, &row_name1);
        if (row_type[0] == 'N') {
            lp->objective_label = string(row_name1);
        }
        else {
			char _type = row_type[0] == 'L' ? '<' : (row_type[0] == 'G' ? '>' : '=');
			Constraint * constraint = new Constraint(_type, 0.0, constr_name_coeff);			
			lp->constraints[row_name1] = constraint;
        }
        return 0;

        case 2:
        row_name2[0] = '\0';
        data2[0] = '0';
        data2[1] = '\0';
        sscanf(line_buf, "%16s %16s %16s %16s %16s", &col_name, &row_name1, &data1, &row_name2, &data2);
        f1 = strtof(data1, NULL);
        f2 = strtof(data2, NULL);

        debug_printf("parser: column: %s %s %s %s %s\n", col_name, row_name1, data1, row_name2, data2);
		
		// It the variable name appears for the first time, then save its name.
		if (lp->variable_names.count(col_name) == 0) {
			lp->variable_names[col_name] = lp->variable_names.size();
		}
		
		// Handle row_name1 / f1.
		// Check if the name refers to the objective function.
		if (strcmp(lp->objective_label.c_str(), row_name1) == 0) {
			lp->objective_name_coeff[col_name] = f1;
		}
		else {
			// The name must correspond to a constraint.
			Constraint * constraint_ptr = lp->constraints[row_name1];
			if (constraint_ptr != NULL) {
				constraint_ptr->name_coeff[col_name] = f1;
			}
			else { printf("#1 NULL pointer!\n"); }
		}

		// Handle row_name2 / f2.
		if (row_name2[0] == '\0') return 0;
		
		// Check if the name refers to the objective function.
		if (strcmp(lp->objective_label.c_str(), row_name2) == 0) {
			lp->objective_name_coeff[col_name] = f2;
		}
		else {
			// The name must correspond to a constraint.
			Constraint * constraint_ptr = lp->constraints[row_name2];
			if (constraint_ptr != NULL) {
				constraint_ptr->name_coeff[col_name] = f2;
			}
			else { printf("#2 NULL pointer!\n"); }
		}
		
        return 0;

        case 3:
        row_name2[0] = '\0';
        data2[0] = '0';
        data2[1] = '\0';
        sscanf(line_buf, "%16s %16s %16s %16s %16s", &rhs_name, &row_name1, &data1, &row_name2, &data2);
        f1 = strtof(data1, NULL);
        f2 = strtof(data2, NULL);
		
		debug_printf("parser: column: %s %s %s %s %s\n", col_name, row_name1, data1, row_name2, data2);
		{
			// Handle row_name1 / f1.
			Constraint * constraint_ptr = lp->constraints[row_name1];
			if (constraint_ptr != NULL) {
				constraint_ptr->rhs = f1;
			}
		
			// Handle row_name2 / f2.
			if (row_name2[0] == '\0') return 0;
			constraint_ptr = lp->constraints[row_name2];
			if (constraint_ptr != NULL) {
				constraint_ptr->rhs = f2;
			}
		}
        return 0;

        case 4:

        sscanf(line_buf, "%16s %16s %16s %16s", &row_type, &row_name1, &col_name, &data1);
		
		debug_printf("parser: column: %s %s %s %s\n", row_type, row_name1, col_name, data1);
		
		f1 = strtof(data1, NULL);
		map<string, double> _constr_name_coeff;
		_constr_name_coeff[col_name] = 1.0;
		char _type = '<';
		switch (row_type[0]) {
			case 'U':
			if (f1 < 0.0) {
				// negative upper bound: we substitute z = -x, and add regular constraint z >= f1
				lp->var_shifts[col_name] = 0.0;
				_type = '>';
				f1 = -f1;
			}
			break;
			case 'L':
			_type = '>';
			if (f1 < 0.0) {
				// negative lower bound: we substitute variable x adding shift: z = x + f1
				lp->var_shifts[col_name] = f1;
				f1 = 0.0;
			}
			break;
		}
		
		// if bound is positive, then a simple constraint is added
		if (f1 > 0.0) {
			Constraint * constraint = new Constraint(_type, f1, _constr_name_coeff);
			lp->constraints[string("#BOUND") + tostr<int>(bounds_count++)] = constraint;
		}
		
        return 0;
    }

    return 0;
}
