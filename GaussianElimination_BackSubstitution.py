import numpy as np

def GaussianElimination(A, y):
    size_mtx = len(A)
    for i in range(size_mtx):
        for j in range(i+1, size_mtx):
            k = A[j][i]/A[i][i]
            A[j] = np.subtract(A[j],np.multiply(k, A[i]))
            y[j] = y[j] - y[i]*k
    return A,y

def BackSubstitution(A, y):
    size_mtx = len(A)
    x = np.zeros(size_mtx)
    x[size_mtx-1] = 1/A[size_mtx-1][size_mtx-1]*y[size_mtx-1]
    p = 0
    for i in range(size_mtx-2,-1,-1):
        for j in range(i+1, size_mtx):
            p += A[i][j]*x[j]
        x[i] = 1/A[i][i] * (y[i] - p)
    return x

def GaussianElimination_BackSubstitution(A, y):
    A, y = GaussianElimination(A, y);
    x = BackSubstitution(A, y);
    return x
    
#A = np.array([[10,5],[2,9]],dtype=float)
#y = np.array([6,3],dtype=float)
A = np.array([[2, 3, -1], [-4, 6, 8], [10,12,14]], dtype=float)
y = np.array([5, 7, 9], dtype=float)

A, y = GaussianElimination(A,y);
x = BackSubstitution(A, y);
print("A", A)
print("y", y)
print("x", x)


k = GaussianElimination_BackSubstitution(A, y)
print("k", k)
