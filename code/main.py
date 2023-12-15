# based on https://github.com/j2kun/simplex-algorithm/blob/main/simplex.py
import heapq
import numpy as np

class SimplexLPSolver:
   def __init__(self, c, A, b):
      self.c = c
      self.A = A
      self.b = b

   def maximize(self):
      self.add_slack_variables()
      self.initialize_simplex_tableau()

      print("Initial tableau:")
      print(np.round(self.tableau))

      while self.can_improve():
         pivot = self.find_pivot_index()
         print("Next pivot index is=%d,%d \n" % pivot)
         self.pivot_about(pivot)
         print("Tableau after pivot:")
         print(self.tableau)
      
      return self.tableau, self.primal_solution(), self.objective_value()

   def add_slack_variables(self):
      self.m, self.n = self.A.shape
      self.A = np.concatenate((self.A, np.eye(self.m)), axis=1)
      self.c = np.concatenate((self.c, np.zeros(self.m + 1))) # +1 since we store the objective function as well

   def initialize_simplex_tableau(self):
      self.tableau = np.concatenate((self.A, self.b.reshape(-1,1)), axis=1)
      self.tableau = np.concatenate((self.tableau, self.c.reshape(1,-1)), axis=0)

   def can_improve(self):
      lastRow = self.tableau[-1]
      return any(map(lambda x : x > 0, lastRow[:-1]))
   
   def find_pivot_index(self):
      # pick minimum positive index of the last row
      column_choices = [(i,x) for (i,x) in enumerate(self.tableau[-1][:-1]) if x > 0]
      column = min(column_choices, key=lambda a: a[1])[0]

      # check if unbounded
      if all(row[column] <= 0 for row in self.tableau):
         raise Exception('Linear program is unbounded.')

      # check for degeneracy: more than one minimizer of the quotient
      quotients = [(i, r[-1] / r[column])
         for i,r in enumerate(self.tableau[:-1]) if r[column] > 0]

      if self.more_than_one_min(quotients):
         raise Exception('Linear program is degenerate.')

      # pick row index minimizing the quotient
      row = min(quotients, key=lambda x: x[1])[0]

      return row, column
   
   def more_than_one_min(self, L):
      if len(L) <= 1:
         return False

      x,y = heapq.nsmallest(2, L, key=lambda x: x[1])
      return x == y
   
   def pivot_about(self, pivot):
      i,j = pivot

      pivotDenom = self.tableau[i][j]
      self.tableau[i] = [x / pivotDenom for x in self.tableau[i]]

      for k,row in enumerate(self.tableau):
         if k != i:
            pivotRowMultiple = [y * self.tableau[k][j] for y in self.tableau[i]]
            self.tableau[k] = [x - y for x,y in zip(self.tableau[k], pivotRowMultiple)]

   def primal_solution(self):
      # the pivot columns denote which variables are used
      columns = np.transpose(self.tableau)
      indices = [j for j, col in enumerate(columns[:-1]) if self.is_pivot_col(col)]
      return [(colIndex, self.variable_value_for_pivot_column(self.tableau, columns[colIndex]))
               for colIndex in indices]
   
   def is_pivot_col(self, col):
      return (len([c for c in col if c == 0]) == len(col) - 1) and sum(col) == 1

   def variable_value_for_pivot_column(self, tableau, column):
      pivotRow = [i for (i, x) in enumerate(column) if x == 1][0]
      return tableau[pivotRow][-1]
   
   def objective_value(self):
      return -self.tableau[-1][-1]

if __name__ == "__main__":
   np.set_printoptions(precision=2)

   A = np.array([[100,   200,    50]
   ,    [5,     12,     11]
   ,    [1,     1,      1]])

   b = np.array([   10**5,  10000,   1000])

   c = np.array([   30,    40,    35])

   solver = SimplexLPSolver(c, A, b)

   tableau, solution, obj = solver.maximize()
   print("Final tableau:")
   print(tableau)
   print("Solution: ", solution)
   print("Objective: ", obj)
   