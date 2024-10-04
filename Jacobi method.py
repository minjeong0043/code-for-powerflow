import numpy as np

def Jacobi(A, y):
    D = np.zeros((len(A), len(A)))
    for i in range(len(A)):
        for j in range(len(A)):
            if i != j:
                D[i][j] = 0
            else:
                D[i][j] = A[i][j]

    M = np.linalg.inv(D) @ (np.subtract(D,A))

    size_mtx = len(A)

    x = np.zeros(size_mtx)
    x_pre = np.zeros(size_mtx)
    e = np.zeros(size_mtx)
    x[0] = 0
    x[1] = 0
    i = 0
    while True:
        x_pre = x.copy()
        x = M @ x + np.linalg.inv(D) @ y
        i += 1
        #print("i", i)
        #print("x", x)
        if np.all(x_pre != 0):
            e[0] = np.abs((x[0] - x_pre[0]) / x_pre[0])
            e[1] = np.abs((x[1] - x_pre[1]) / x_pre[1])
            if np.all(e < 10 ** (-4)):
                print("e = ", e)
                print("i = ", i)
                break
    return x

A = np.array([[10,5],[2,9]],dtype=float)
y = np.array([6,3],dtype=float)

x = Jacobi(A, y)