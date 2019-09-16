import sys, os
import subprocess
from pathlib import Path
from xml.dom import minidom
from utils import *

SCRIPT_PATH = os.path.dirname(Path(__file__).absolute())
DATA_DIR = os.path.join(SCRIPT_PATH, "..", "data")

EXEC_PATH = os.path.join(SCRIPT_PATH, "..", "simplex")
CPLEX_PATH = "/Users/maciek/Applications/IBM/ILOG/CPLEX_Studio127/cplex/bin/x86-64_osx/cplex"

def extract_sol_from_xml(path):
    sol = dict()
    xmldoc = minidom.parse(path)
    itemlist = xmldoc.getElementsByTagName("variable")
    for item in itemlist:
        sol[item.attributes["name"].value] = float(item.attributes["value"].value)

    objective = float(xmldoc.getElementsByTagName("header")[0].attributes["objectiveValue"].value)
    
    return sol, objective

def test_compare_with_cplex():
    success_count = 0
    failure_count = 0
    test_cases = []
    for filename in get_files_from_dir(DATA_DIR, ext = ".lp"):
        call_list = [EXEC_PATH, filename]
        subprocess.call(call_list)
        sol_path = os.path.join(SCRIPT_PATH, "sol_test.xml")
        sol1, objective1 = extract_sol_from_xml(sol_path)
        
        sol_path = os.path.join(SCRIPT_PATH, "sol_cplex.xml")
        call_list = [CPLEX_PATH, "-c", "read", filename, "opt", "write", sol_path, "sol", "y"]
        subprocess.call(call_list)
        sol2, objective2 = extract_sol_from_xml(sol_path)
        
        if set(sol1.keys()) != set(sol2.keys()):
            print("Problems have different sets of variables!")
        else:
            equal_sols = True    
            for var, val in sol1.items():
                if abs(sol2[var] - val) > 1e-3:
                    print("Different solutions for Test and Cplex.")
                    equal_sols = False
                    break
        
            if equal_sols:
                print("Solutions of Test and Cplex are the same.")
        
        print("Objective value: Test:", objective1, "Cplex:", objective2)
        
        if abs(objective1 - objective2) < 1e-3:
            success_count += 1
        else:
            failure_count += 1
            
        test_cases.append((os.path.basename(filename), objective1, objective2))
        
    print("\nSuccess Count:", success_count, "\nFailure Count:", failure_count)
    print(test_cases)
        
test_compare_with_cplex()

