# note: for "rhs" values less than 4 the iteration number is much lower
from sys import argv

n = int(argv[1])
rhs = 5

with open("./data/klee-minty-n" + str(n) + ".lp", "w") as f:
    print("Maximize", end='\n', file=f)
    for i in range(1,n+1):
        print(str(2**(n-i)) + "x" + str(i) + ("+" if i < n else ""), end='', file=f)
    print("\nSubject To", end='\n', file=f)
    print("c1:  x1 <=" + str(rhs), end='\n', file=f) 
    for i in range(2,n+1):
        print("c" + str(i) + ": ", end='', file=f)
        for j in range(1, i):
            print(str(2**(i-j+1)) + "x" + str(j) + "+", end='', file=f)
        print("x" + str(i) + "<=" + str(rhs**i), end='\n', file=f)
    print("End\n", file=f)

 