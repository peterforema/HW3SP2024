import copy

def symmetric(matrix):
    n = len(matrix)
    for i in range(n):
        for j in range(i + 1, n):
            if matrix[i][j] != matrix[j][i]:
                return False
    return True

def positive_definite(matrix):
    n = len(matrix)
    for i in range(n):
        det = matrix[i][i]
        for j in range(i):
            det -= matrix[i][j] ** 2
        if det <= 0:
            return False
    return True

def cholesky_method(A):
    n = len(A)
    L = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            if i == j:
                sum_val = sum(L[i][k] ** 2 for k in range(j))
                L[i][j] = (A[i][j] - sum_val) ** 0.5
            else:
                sum_val = sum(L[i][k] * L[j][k] for k in range(j))
                L[i][j] = (A[i][j] - sum_val) / L[j][j]
    b = [r.pop() for r in A]
    y = forward_substitution(L, b)
    x = backward_substitution(transpose(L), y)
    return x

def forward_substitution(L, b):
    # some code adapted from ChatGPT
    '''function takes 2 arguments : lower matrix L and vector b
    empty vector is y and it is the smame length as L and b
    '''
    n = len(L)
    y = [0] * n
    for i in range(n):
        y[i] = b[i]
        for j in range(i):
            y[i] -= L[i][j] * y[j]
        y[i] /= L[i][i]
    return y

def backward_substitution(LT, y):
    # code adapted from ChatGPT
    '''this function takes 2 arguments, upper matrix LT and vector y
    empty vector is x and it is the same length as LT and y
    This function iterates in reverse order of the previous function'''
    n = len(LT)
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = y[i]
        for j in range(i + 1, n):
            x[i] -= LT[i][j] * x[j]
        x[i] /= LT[i][i]
    return x

def transpose(matrix):
    '''transform horizontal matrix A into vertical matrix AT
    this is can be helpful for factorization'''
    return [list(row) for row in zip(*matrix)]



def LUFactorization(A):
    '''Copied from Prof Smays DoolittleMethod.py Example'''
    """
    This is the Lower-Upper factorization part of Doolittle's method.  The factorizaiton follows the work in
    Kreyszig section 20.2.  Note: L is the lower triangular matrix with 1's on the diagonal.  U is the upper traingular matrix.
    :param A: a nxn matrix
    :return: a tuple with (L, U)
    """
    n = len(A)
    # Step 1
    U = [([0 for c in range(n)] if not r == 0 else [a for a in A[0]]) for r in range(n)]
    L = [[(1 if c==r else (A[r][0]/U[0][0] if c==0 else 0)) for c in range(n)] for r in range(n)]

    #step 2
    for j in range(1,n):  # j is row index
        #(a)
        for k in range(j,n):  # always j >= 1 (i.e., second row and higher)
            U[j][k]=A[j][k]  # k is column index and scans from column j to n-1
            for s in range(j):  #  s is column index for L and row index for U
                U[j][k] -= L[j][s]*U[s][k]
            #(b)
            for i in range(k+1, n):
                sig=0
                for s in range(k):
                    sig+=L[i][s]*U[s][k]
                L[i][k]=(1/(U[k][k]))*(A[i][k]-sig)
    return (L,U)

def BackSolve(A,b,UT=True):
    '''Copied from Prof Smays DoolittleMethod.py Example'''

    """
    This is a backsolving algorithm for a matrix and b vector where A is triangular
    :param A: A triangularized matrix (Upper or Lower)
    :param b: the right hand side of a matrix equation Ax=b
    :param UT: boolean of upper triangular (True) or lower triangular (False)
    :return: the solution vector x, from Ax=b
    """
    nRows=len(b)
    nCols=nRows
    x=[0]*nRows
    if UT:
        for nR in range(nRows-1,-1,-1):
            s=0
            for nC in range(nR+1,nRows):
                s+=A[nR][nC]*x[nC]
            x[nR]=1/A[nR][nR]*(b[nR]-s)
    else:
        for nR in range(nRows):
            s=0
            for nC in range(nR):
                s+=A[nR][nC]*x[nC]
            x[nR]=1/A[nR][nR]*(b[nR]-s)
    B = GS.checkMatrixSoln(A, x, False)
    return x

def Doolittle(Aaug):
    '''Copied from Prof Smays DoolittleMethod.py Example'''

    """

    :param Aaug: the augmented matrix
    :return: the solution vector x
    """
    A,b=GS.separateAugmented(Aaug)
    L,U=LUFactorization(A)
    B=GS.matrixMult(L,U)
    y=BackSolve(L,b, UT=False)
    x=BackSolve(U,y, UT=True)
    return x

def solve_matrix_equation(A):
    '''trying Cholesky's method before using Doolittle if necessary'''
    try:
        x = cholesky_method(A)
        method_used = "Cholesky"
    except ValueError:
        x = Doolittle(A)
        method_used = "Doolittle"
    return x, method_used


def main():
    # Problem 1
    A1 = [
        [1, -1, 3, 2, 15],
        [-1, 5, -5, -2, -35],
        [3, -5, 19, 3, 94],
        [2, -2, 3, 21, 1]
    ]
    x1, method_used1 = solve_matrix_equation(A1)
    print(f"Problem 1 - Using {method_used1} method: {x1}")

    # Problem 2
    A2 = [
        [4, 2, 4, 0, 20],
        [2, 2, 3, 2, 36],
        [4, 3, 6, 3, 60],
        [0, 2, 3, 9, 122]

    ]
    x2, method_used2 = solve_matrix_equation(A2)
    print(f"Problem 2 - Using {method_used2} method: {x2}")

if __name__ == "__main__":
    main()
