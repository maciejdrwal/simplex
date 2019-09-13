/ value = 4.5
/ x* = [0, 0.5, 0, 2.5, 1.5]
/ NOTE: if we start from basis (0,1,2) we end up with infeasible solution

Minimize
 obj: x1 + x2 + x3 + x4 + x5
 
Subject To
 C0: 3 x1 + 2 x2 + x3           = 1
 C1: 5 x1 +   x2 + x3 + x4      = 3
 C2: 2 x1 + 5 x2 + x3      + x5 = 4
End
