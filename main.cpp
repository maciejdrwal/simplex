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
    
    if (argc < 2)
    {
        std::cout << "You must provide a file name." << std::endl;
        return 1;
    }

    simplex::FileFormat file_format;

    if (strcmp(strrchr(argv[1], '.'), ".lp") == 0) 
    {
        file_format = simplex::FileFormat::LP_FILE;
    }
    else if (strcmp(strrchr(argv[1], '.'), ".mps") == 0) 
    {
        file_format = simplex::FileFormat::MPS_FILE;
    }
    else 
    {
        std::cout << "Unrecognized file extension (must be one of: .lp, .mps)." << std::endl;
        return 1;
    }
    
    std::ifstream fin;
    fin.open(argv[1]);
    
    if (!fin) 
    {
        std::cout << "Invalid input file name." << std::endl;
        return 1;
    }

    // Read whole input file instead of separate lines
    fin.seekg(0, std::ios::end);
    const size_t size = fin.tellg();
    std::string buffer(size, ' ');
    fin.seekg(0);
    fin.read(&buffer[0], size);
    fin.close();

    // Alternative method:
    //std::stringstream buffer;
    //buffer << fin.rdbuf();

    std::cout << "Reading file done.\n" << std::endl;

    simplex::LinearProgram lp;
    simplex::Presolve presolve(lp);

    try 
    {
        if (simplex::run_parser_lp(buffer, lp))
        {
            presolve.run();
            lp.solve();
            lp.write("sol_test.xml");
        }
    }
    catch (const std::string & msg) 
    {
        std::cout << "Exception: " << msg << std::endl;
    }
    
    std::cout << "Exiting program." << std::endl;
    
    return 0;
}
