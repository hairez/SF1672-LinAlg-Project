from scipy import optimize as sp



A = [   [100,   200,    50]
    ,   [5,     12,     11]
    ,   [1,     1,      1]]

b = [   10**5,  10000,   1000]


c = [   -30,    -40,    -35]


res = sp.linprog(c, A_ub = A, b_ub = b, method  = 'simplex')

print(res.fun)
print(res.x)
