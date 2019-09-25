/ x* = [0, 1, 3, 0, 2, 0, 0]
/ value = 2

Minimize
 obj: 0 x1 + 2 x2 + 0 x3 + x4 + 0 x5 + 0 x6 + 5 x7
 
Subject To
 C0: x1 + x2 + x3 + x4                = 4
 C1: x1                + x5           = 2
 C2:           x3           + x6      = 3
 C3:    3 x2 + x3                + x7 = 6
Bounds
 1 <= x2 <= 10
 0 <= x1
 2 <= x5
 2.5 <= x3 <= 4
 x4 <= 0.5
End
