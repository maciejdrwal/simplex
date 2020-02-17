// Copyright (C) 2019 Maciej Drwal
// 
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
// 

#include <iostream>
#include <fstream>

#include "parser.h"
#include "presolve.h"
#include "simplex.h"

int main(int argc, char ** argv) 
{
    std::cout.precision(4);
    std::cout.setf(std::ios_base::showpoint);
    
    if (argc < 2) {
        std::cout << "You must provide a file name." << std::endl;
        return 1;
    }

    if (strcmp(strrchr(argv[1], '.'), ".lp") == 0) {
        file_format = LP_FILE;
    }
    else if (strcmp(strrchr(argv[1], '.'), ".mps") == 0) {
        file_format = MPS_FILE;
    }
    else {
        std::cout << "Unrecognized file extension (must be one of: .lp, .mps)." << std::endl;
        return 1;
    }
    
    std::ifstream fin;
    fin.open(argv[1]);
    
    if (!fin) {
        std::cout << "Invalid input file name." << std::endl;
        return 1;
    }


    // Read whole input file instead of separate lines
    fin.seekg(0, std::ios::end);
    size_t size = fin.tellg();
    std::string buffer(size, ' ');
    fin.seekg(0);
    fin.read(&buffer[0], size);
    fin.close();

    // Alternative method:
    //std::stringstream buffer;
    //buffer << fin.rdbuf();

    std::cout << "Reading file done.\n" << std::endl;

    LinearProgram lp;
    Presolve presolve(lp);

    // Code used to feed old hand-written parsers (non-Spirit)
    // string line;
    // while (fin) {        
    //     getline(fin, line);
    //     line += "\n";
        
    //     if (file_format == LP_FILE) {            
    //         if (parse_input_line_lp(line.c_str(), &lp) != 0) {
    //             fin.close();
    //             std::cout << "Aborting." << std::endl;
    //             return 1;
    //         }
    //     }
    //     else if (file_format == MPS_FILE) {
    //         if (parse_input_line_mps(line.c_str(), &lp) != 0) {
    //             fin.close();
    //             std::cout << "Aborting." << std::endl;
    //             return 1;
    //         }
    //     }
    // }

    try {
        if (run_parser_lp(buffer, lp)) {
            presolve.run();
            lp.solve();
            lp.write("sol_test.xml");
        }
    }
    catch (const std::string& msg) {
        std::cout << "Exception: " << msg << std::endl;
    }
    
    std::cout << "Exiting program." << std::endl;
    
    return 0;
}
