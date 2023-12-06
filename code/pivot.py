eps = 1e-8
inf = 2e9

class LPSolver:
    
    def pivot(rr,s):
        a = D[r][:]
        inv = 1/a[s]
        for i in range(m+2):
            if i != r and abs(D[i][s] > eps):
                b = D[i][:]
                inv2 = b[s]*inv

        for j in range(n+2):
            if j != s: D[r][j] *= inv
        for i in range(m+2):
            if i != r: D[i][s] *= -inv
        D[r][s] = inv
        B[r],N[s] = N[s],B[r]


