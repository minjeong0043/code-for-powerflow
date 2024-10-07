import numpy as np
from sympy import Matrix, symbols


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

#x1, x2 = symbols('x1, x2')
#Newton_Raphson(Matrix([x1+x2, x1*x2]), Matrix([x1, x2]), Matrix([4, 9]), Matrix([15, 50]))

x = symbols('x')
Newton_Raphson(Matrix([x**2]), Matrix([x]), Matrix([1]), Matrix([9]))
