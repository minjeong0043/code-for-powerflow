import numpy as np
from sympy import Matrix, symbols

def GaussianElimination(A, y):
    # print("A" ,A)
    # print("Y", y)
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

def Jacobi_GaussSeidel(A, y):
    D = np.zeros((len(A), len(A)))
    for i in range(len(A)):
        for j in range(len(A)):
            if i != j: # Jacobi
            #if i < j: # Gauss-Seidel
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

        if np.all(x_pre != 0):
            e[0] = np.abs((x[0] - x_pre[0]) / x_pre[0])
            e[1] = np.abs((x[1] - x_pre[1]) / x_pre[1])
            if np.all(e < 10 ** (-4)):
                print("e = ", e)
                print("i = ", i)
                break
    return x

def Newton_Raphson(f, x, x_init, y):
    # Jacobian Mtx
    size_mtx = len(f)
    J = Matrix(size_mtx, size_mtx, lambda i, j: f[i].diff(x[j]))
    #J_inv = J.inv()

    x_val = x_init.copy()
    e = np.zeros(size_mtx)
    i  = 0

    while True:
        x_pre = x_val.copy()
        x_val_np = np.array(x_val, dtype=float)
        print(x_val_np)

        J_eval = J.subs({x[i]: x_pre[i] for i in range(size_mtx)})
        J_inv = J_eval.inv()
        f_eval = Matrix([func.subs({x[i]: x_pre[i] for i in range(size_mtx)}) for func in f])
        #print("f_eval : ",f_eval)
        #print("J_inv : ", J_inv)
        x_val = x_pre + J_inv @ (y - f_eval)
        #x_val_np = np.array(x_val, dtype=float)
        #print(x_val_np)


        e = np.zeros(size_mtx)
        for k in range(size_mtx):
            e[k] = np.abs((x_val[k] - x_pre[k]) / x_pre[k])
        i += 1
        x_pre = x_val
        if np.all(x_pre != 0):
            if np.all(np.array(e, dtype=float) < 10 ** (-4)):
                print("i = ", i)
                print("e = ", np.array(e, dtype = float))
                print("x = ", np.array(x_val, dtype=float))
                break

        if i > 1000:
            print("Didn't converge until interation num 1000!!!!")
            break

# def Newton_Raphson_4phase(f, x, x_pre, y):
def Newton_Raphson_4phase(f, x, x_init, y):
    # Jacobian Mtx
    size_mtx = len(f)
    J = Matrix(size_mtx, size_mtx, lambda i, j: f[i].diff(x[j]))
    #J_inv = J.inv()
    i = 0

    x_val = x_init.copy()

    del_y = 0
    del_x = 0

    while True:
        # Phase 1
        del_y = y - f.subs({x[i]: x[i] for i in range(size_mtx)})
        del_y = del_y.subs({x[i]: x_val[i] for i in range(size_mtx)})
        del_y = np.array(del_y, dtype=float)

        # Phase 2
        J = Matrix(size_mtx, size_mtx, lambda i, j:f[i].diff(x[j]))
        J_val = J.subs({x[i]: x_val[i] for i in range(size_mtx)})
        J_val = np.array(J_val, dtype=float)

        # Phase 3
        del_x = GaussianElimination_BackSubstitution(J_val, del_y)
        del_x = Matrix(del_x)

        # Phase 4
        x_pre = np.array(x_val, dtype=float)
        # print("x_pre : ", x_pre)
        # print("del_x : ", del_x)
        x_val = x_pre + del_x
        # print("x_val : ", x_val)

        e = np.zeros(size_mtx)
        for k in range(size_mtx):
            e[k] = np.abs((x_val[k] - x_pre[k]) / x_pre[k])
        i += 1
        if np.all(x_pre != 0):
            if np.all(np.array(e, dtype=float) < 10 ** (-4)):
                print("i = ", i)
                print("e = ", np.array(e, dtype = float))
                print("x = ", np.array(x_val, dtype=float))
                break

        if i > 1000:
            print("Didn't converge until interation num 1000!!!!")
            break