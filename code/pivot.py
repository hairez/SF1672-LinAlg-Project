#Code is based off KACTL's implementation of the Simplex algorithm
#which has its source from the Stanford Notebook

#https://github.com/kth-competitive-programming/kactl/blob/81d547a917d05fd482ba2e55a3d1fb1e444af919/content/numerical/Simplex.h

eps = 1e-8
inf = 2e9



class LPSolver:
    def __init__(self, A, b, c):

        self.m = len(b)
        self.n = len(c)
        self.N = [0]*(self.n+1) #array med längd n+1
        self.B = [0]*(self.m)
        self.D = [[0]*(self.n+2) for _ in range(self.m+2)]

        for i in range(self.m):
            for j in range(self.n):
                self.D[i][j] = A[i][j]
        
        for i in range(self.m):
            self.B[i]=self.n+i
            self.D[i][self.n] = -1
            self.D[i][self.n+1] = b[i]
        
        for j in range(self.n):
            self.N[j] = j
            self.D[self.m][j] = - c[j]
        
        self.N[self.n] = -1
        self.D[self.m+1][self.n] = 1
    
    def pivot(self, r, s):
        a = self.D[r]
        inv = 1/a[s]
        for i in range(self.m+2):
            if i != r and abs(self.D[i][s]) > eps:
                b = self.D[i]
                inv2 = b[s]*inv

                for j in range(self.n+2):
                    b[j] -= a[j] * inv2
                b[s] = a[s] * inv2


        for j in range(self.n+2):
            if j != s: self.D[r][j] *= inv
        for i in range(self.m+2):
            if i != r: self.D[i][s] *= -inv
        self.D[r][s] = inv
        self.B[r], self.N[s] = self.N[s], self.B[r]
    
    def ltj(self, X, N, s, j):
        if s == -1 or (X[j], N[j]) < (X[s], N[s]):
            return j
        return s
    
    def simplex(self, phase):
        x = self.m + phase -1
        while 1:
            s = -1
            for j in range(self.n+1):
                if self.N[j] != -phase:
                    s = self.ltj(self.D[x], self.N, s, j)
            if self.D[x][s] >= -eps:return 1

            r = -1
            for i in range(self.m):
                if self.D[i][s] <= eps:
                    continue
                if r == -1 or (self.D[i][self.n+1]/self.D[i][s], self.B[i]) < (self.D[r][self.n+1]/self.D[r][s], self.B[r]):
                    r = i

            if r ==-1:
                return 0
            
            self.pivot(r,s)
    
    def solve(self):
        r = 0
        for i in range(1,self.m):
            if self.D[i][self.n+1] < self.D[r][self.n+1]:
                r = i
        
        if self.D[r][self.n+1] < -eps:
            self.pivot(r,self.n)
            if (not self.simplex(2)) or self.D[self.m+1][self.n+1] < -eps:
                return -inf
            for i in range(self.m):
                if self.B[i] == -1:
                    s = 0
                    for j in range(1,self.n+1):
                        s = self.ltj(self.D[i], self.N, s, j)
                    self.pivot(i, s)

        ok = self.simplex(1)

        x = self.x = [0]*(self.n)
        for i in range(self.m):
            if self.B[i] < self.n:
                x[self.B[i]] = self.D[i][self.n+1]
        return self.D[self.m][self.n+1] if ok else inf
                



