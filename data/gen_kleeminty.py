# note: for "rhs" values less than 4 the iteration number is much lower
from sys import argv
n = int(argv[1])
rhs = 5
f = open("./data/klee-minty-n" + str(n) + ".lp", "wb")
print >>f, "Maximize"
for i in xrange(1,n+1):
    print >>f, 2**(n-i), "x" + str(i), "+" if i < n else "",
print >>f, "\nSubject To"
print >>f, "c1:  x1 <=", rhs
for i in xrange(2,n+1):
    print >>f, "c" + str(i) + ": ",
    for j in xrange(1, i):
        print >>f, 2**(i-j+1), "x" + str(j), "+",
    print >>f, "x" + str(i), "<=", rhs**i
print >>f, "End\n"

 