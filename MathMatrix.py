'''
Noé de Lima Bezerra
30/12/2016 - 2h53

Testes de Funções para utilização posterior
Operações com matrizes
'''
import sys


#Verify if a matrix is square
def is_square(A):
    if not type(A) == list:
        print('Is not a matrix!')
        return False
    else:
        size = len(A)
        for i in range(size):
            if not type(A[i]) == list:
                print("Dimension error: doesn't have a valid column")
                return False
            elif not len(A[i]) == size:
                print('Dimension error: not a square matrix!')
                return False
            else:
                for j in range(size):
                    if type(A[i][j]) == list:
                        print('Dimension error: multidimensional element in matrix!')
                        return False
    return True


#Criate a diagonal matrix with the elements of a vector
def diag(V):
    D = []
    for i in range(len(V)):
        D.append([])
        if type(V[i]) == list:
            print('Vector with multidimensional value!')
            return 0
        for j in range(len(V)):
            if i == j:
                D[i].append(V[i])
            else:
                D[i].append(0)
    return D


#Calculate determinant of a matrix A
def det(A):
    #Verify all conditions
    key = is_square(A)
    if not key:
        print('Determinant error: the square format verification return false')
        return None
    out = 0  #Initialize det as 0 value
    size = len(A)
    if size == 1:
        return A[0][0]
    for j in range(size):  #Calculating Determinant
        M = []
        a = A[0][j]
        signal = ((-1) ** j)
        if not a == 0:
            for line in range(1, size):
                M.append([])
                for row in range(size):
                    if not row == j:
                        M[line-1].append(A[line][row])                        
            out += signal*a*det(M)
            #print('The current value of det for', size, '-order array is', out, ', and will add', signal*a, 'to the next det')
    #print('The final value of det is', out)
    return out


#Multiply two matrix M and N
def multiply(M, N):
    m_line = len(M)
    m_row = len(M[0])
    n_line = len(N)
    n_row = len(N[0])

    #Verify all conditions
    if not m_row == n_line:
        print('Dimension error!')
        return None
    for i in range(m_line):
        if not len(M[i]) == m_row:
            print("The length of M isn't uniform")
            return None
        for j in range(m_row):
            if type(M[i][j]) == list:
                print('Multidimensional element in matrix M')
                return None
    for i in range(n_line):
        if not len(N[i]) == n_row:
            print("The length of N isn't uniform")
            return None
        for j in range(n_row):
            if type(N[i][j]) == list:
                print('Multidimensional element in matrix N')
                return None
            
    P = []#Create P as product matrix to calculate
    for i in range(m_line):
        P.append([])
        for j in range(n_row):
            element = 0
            for k in range(m_row):
                element += M[i][k]*N[k][j]
            P[i].append(element)
    return P


#Calculate the inverse of a matrix
def inverse(A):
    key = is_square(A)
    #print(key)
    if not key:  #Verify all conditions
        print('Inverse error: not a valid square matrix')
        return None
    elif det(A) == 0:
        print("Singular matrix: can't be inverted")
        return(None)
    
    size = len(A)  #The size of matrix
    I = []  #Define a matrix I as A, to transform to a ident matrix
    for i in range(size):
        I.append([])
        for j in range(size):
            I[i].append(A[i][j])
    V = []  #Create a ident matrix to transforma to a inverse matrix
    for i in range(size):
        V.append(1)
    R = diag(V)
    #Don't let be zero in divisor addicting another nonzero line
    for i in range(size):
        k = 0
        while I[i][i] == 0:
            k += 1
            if not k == i:
                for row in range(size):
                    I[i][row] += I[k][row]
                    R[i][row] += R[k][row]

        a = I[i][i]
        for row in range(size):  #Divide the line for main diagonal element
            I[i][row] = I[i][row]/a
            R[i][row] = R[i][row]/a
        for line in range(size):  #Do the anothers columns elements equal zero
            if not line == i:
                a = I[line][i]
                for row in range(size):
                    I[line][row] -= a*I[i][row]
                    R[line][row] -= a*R[i][row]
    return R


def transpose(A):
    a_line = len(A)
    a_row = len(A[0])
    #Verify conditions
    for i in range(a_line):
        if not len(A[i]) == a_row:
            print("The length of A isn't uniform")
            return None
        for j in range(a_row):
            if type(A[i][j]) == list:
                print('Multidimension element in matrix A')
                return None
    #Calculating the transposed matrix
    T = []
    for i in range(a_row):
        T.append([])
        for j in range(a_line):
            T[i].append(A[j][i])
    return T


#Function to print Matrix
def printmatrix(A, fh=sys.stdout):
    for i in range(len(A)):
        print('Line', i+1, ':', end='', file=fh)
        for j in range(len(A[i])):
            print('\t', A[i][j], end='', file=fh)
        print('', file=fh)
