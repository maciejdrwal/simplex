// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#include "Parser.h"
#include "Presolve.h"
#include "Simplex.h"

#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{
    std::cout.precision(4);
    std::cout.setf(std::ios_base::showpoint);

    if (argc < 2)
    {
        std::cout << "You must provide a file name." << std::endl;
        return 1;
    }

    const std::string_view filename(argv[1]);
    if (filename.substr(filename.find_last_of('.') + 1) != "lp")
    {
        std::cout << "Unrecognized file extension (must be .lp)." << std::endl;
        return 1;
    }

    std::ifstream fin;
    fin.open(argv[1]);

    if (!fin)
    {
        std::cout << "Invalid input file name." << std::endl;
        return 1;
    }

    // Read whole input file
    fin.seekg(0, std::ios::end);
    const auto size = fin.tellg();
    std::string buffer(size, ' ');
    fin.seekg(0);
    fin.read(buffer.data(), size);
    fin.close();

    std::cout << "Reading file done.\n";

    try
    {
        simplex::LinearProgram lp;
        run_parser_lp(buffer, lp);
        buffer.clear();
        simplex::Presolve presolve(lp);
        presolve.run();
        lp.solve();
        lp.write("sol_test.xml");
    }
    catch (const std::string &msg)
    {
        std::cout << "Exception: " << msg << std::endl;
    }

    std::cout << "Exiting program." << std::endl;

    return 0;
}
