import sys, os
import subprocess
import platform
import cplex
from pathlib import Path
from xml.dom import minidom
from utils import *

SCRIPT_PATH = os.path.dirname(Path(__file__).absolute())
DATA_DIR = os.path.join(SCRIPT_PATH, "..", "data", "test_set_1")
EXEC_PATH_UNIX = os.path.join(SCRIPT_PATH, "..", "simplex")
EXEC_PATH_WIN = os.path.join(SCRIPT_PATH, "..\\build\\Debug",  "simplex.exe")
EXEC_PATH = EXEC_PATH_WIN if platform.system()=="Windows" else EXEC_PATH_UNIX

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
    passed_tests = []
    failed_tests = []
    for filename in get_files_from_dir(DATA_DIR, ext = ".lp"):
        print(">>> ", filename)
        call_list = [EXEC_PATH, filename]
        print(call_list)
        subprocess.call(call_list)
        sol_path = os.path.join(SCRIPT_PATH, "sol_test.xml")
        sol1, objective1 = extract_sol_from_xml(sol_path)
        
        sol_path = os.path.join(SCRIPT_PATH, "sol_cplex.xml")
        
        try:
            c = cplex.Cplex()
            c.read(filename)
            c.solve()
            c.solution.write(sol_path)
            sol2, objective2 = extract_sol_from_xml(sol_path)
        except:
            print("Unable to solve via CPLEX")
            continue
            
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

        test_case = {"name": os.path.basename(filename), "value": objective1, "expected": objective2}
        
        if abs(objective1 - objective2) < 1e-3:
            success_count += 1
            passed_tests.append(test_case)
        else:
            failure_count += 1
            failed_tests.append(test_case)
        
    print("\nSuccess Count:", success_count, "\nFailure Count:", failure_count)
    print("Passed Tests:", passed_tests)
    print("\nFailed Tests:", failed_tests)
        
test_compare_with_cplex()

