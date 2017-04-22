"""
Struct to manipulate numerical matrix objects
Author: Noe de Lima Bezerra
Company: UFF - Universidade Federal Fluminense
Author Email: noe_lima@id.uff.br
Created at 10/03/2017
"""

"""
Define matrix as a class and setting 2 dimension
"""
import sys

class matrix:
    def __init__(self, a, b=None):
        if not b:
            b = a
        self.matrix = []
        for i in range(a):
            self.matrix.append([])
            for j in range(b):
                self.matrix[i].append(0.0)
        self.line = a
        self.row = b

    #A set of matrix functions to general usage

    def printmatrix(self, fh=sys.stdout): #Function to print a matrix in a standard output
        for i in range(self.line):
            print('Line', i+1, ':', end='', file=fh)
            for j in range(self.row):
                print('\t', self.matrix[i][j], end='', file=fh)
            print('', file=fh)
        
    def det(self): #Define determinant function to square matrix
        if not self.line == self.row: #Void return to a non square matrix
            return None
        out = 0  #Initialize det as 0 value
        size = self.line
        if size == 1:
            return self.matrix[0][0]
        for j in range(size):  #Calculating Determinant
            M = matrix(size-1)
            a = self.matrix[0][j]
            signal = ((-1) ** j)
            if not a == 0:
                for line in range(1, size):
                    for row in range(size):
                        if row < j:
                            M.matrix[line-1][row] = self.matrix[line][row]
                        elif row > j:
                            M.matrix[line-1][row-1] = self.matrix[line][row]
                out += signal*a*M.det()
        return out

    def inverse(self): #Function to determinate the inverse of a square matrix
        if not self.line == self.row:
            print("Isn't a square matrix")
            return None #Test if is square
        elif self.det() == 0:
            print("Singular matrix: can't be inverted")
            return None #Test to singularity
        size = self.line  #The size of matrix
        I = matrix(size)  #Define a matrix I as self, to transform to a identity matrix
        for i in range(size):
            for j in range(size):
                I.matrix[i][j] = self.matrix[i][j]
        R = matrix(size)  #Create a ident matrix to transforma to a inverse matrix
        for i in range(size):
            R.matrix[i][i] = 1
        for i in range(size): #Don't let be zero in divisor addicting another nonzero line
            k = 0
            while I.matrix[i][i] == 0:
                k += 1
                if not k == i:
                    for row in range(size):
                        I.matrix[i][row] += I.matrix[k][row]
                        R.matrix[i][row] += R.matrix[k][row]

            a = I.matrix[i][i]
            for row in range(size):  #Divide the line for main diagonal element
                I.matrix[i][row] = I.matrix[i][row]/a
                R.matrix[i][row] = R.matrix[i][row]/a
            for line in range(size):  #Do the anothers columns elements equal zero
                if not line == i:
                    a = I.matrix[line][i]
                    for row in range(size):
                        I.matrix[line][row] -= a*I.matrix[i][row]
                        R.matrix[line][row] -= a*R.matrix[i][row]
        return R
    
    def transpose(self): #Function to transpose a matrix
        T = matrix(self.row, self.line)
        for i in range(self.row):
            for j in range(self.line):
                T.matrix[i][j] = self.matrix[j][i]
        return T
    
    def multiply(self, N): #Function to multiply matrix
        if not type(N) == matrix:
            return None #Verify if N is a matrix too
        if not self.row == N.line: #Verify dimensions
            print('Dimension error!')
            return None
        P = matrix(self.line, N.row) #Create P as product matrix to calculate
        for k in range(self.row): #Using kij sequence to optimize cache memory access
            for i in range(self.line):
                r = self.matrix[i][k]
                for j in range(N.row):
                    P.matrix[i][j] += r*N.matrix[k][j]
        return P

    def multiplysc(self, N): #Multiply matrix by a scalar
        M = matrix(self.line, self.row)
        for i in range(self.line):
            for j in range(self.row):
                M.matrix[i][j] = N*self.matrix[i][j]
        return M

    def matrixsum(self, M): #Calculate de sum of two matrix with the same dimension
        if not type(M) == matrix and self.line == M.line and self.row == M.row:
            return None
        S = matrix(self.line, self.row)
        for i in range(self.line):
            for j in range(self.row):
                S.matrix[i][j] = self.matrix[i][j] + M.matrix[i][j]
        return S
