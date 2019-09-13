// Copyright (C) 2017-2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#include <iostream>
#include <fstream>

#include "parser.h"
#include "simplex.h"

using namespace std;

int main(int argc, char ** argv) 
{
    if (argc < 2) {
        cout << "You must provide a file name." << endl;
        return 1;
    }
    
    if (strcmp(strrchr(argv[1], '.'), ".lp") == 0) {
        file_format = LP_FILE;
    }
    else if (strcmp(strrchr(argv[1], '.'), ".mps") == 0) {
        file_format = MPS_FILE;
    }
    else {
        cout << "Unrecognized file extension (must be one of: .lp, .mps)." << endl;
        return 1;
    }
    
    LinearProgram lp;
    
    ifstream fin;
    fin.open(argv[1]);
    
    if (!fin) {
        cout << "Invalid input file name." << endl;
        return 1;
    }
        
    string line;
    
    while (fin) {        
        getline(fin, line);
        line += "\n";
        
        if (file_format == LP_FILE) {            
            if (parse_input_line_lp(line.c_str(), &lp) != 0) {
                fin.close();
                cout << "Aborting." << endl;
                return 1;
            }
		}
        else if (file_format == MPS_FILE) {
            if (parse_input_line_mps(line.c_str(), &lp) != 0) {
                fin.close();
                cout << "Aborting." << endl;
                return 1;
            }
        }
    }

    fin.close();
    cout << "Reading file done.\n" << endl;
	
	lp.solve();
    lp.write("sol_test.xml");
	cout << "Exiting program." << endl;
	
    return 0;
}
