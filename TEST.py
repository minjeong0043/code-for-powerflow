import numpy as np
from sympy import Matrix, symbols
from Func import GaussianElimination_BackSubstitution

x1, x2 = symbols('x1 x2')
x = Matrix([x1, x2])
f = Matrix([x1+x2, x1*x2])
# x = np.array([x1, x2])
# f = np.array([x1+x2, x1*x2])
size_mtx = len(f)

# initial values
y = Matrix([15, 50])
x_val = Matrix([4, 9])
del_y = 0
del_x = 0

# phase 1
# y_pre = y.copy()
del_y = y - f.subs({x1:x[0], x2:x[1]})
del_y = del_y.subs({x[i]:x_val[i] for i in range(size_mtx)})
del_y = np.array(del_y, dtype=float)
# print(del_y)

# phase 2
J = Matrix(size_mtx, size_mtx, lambda i, j: f[i].diff(x[j]))
J_val = J.subs({x[i]: x_val[i] for i in range(size_mtx)})
J_val = np.array(J_val, dtype=float)
# print(J_val)

# phase 3
del_x = GaussianElimination_BackSubstitution(J_val, del_y)
del_x = Matrix(del_x)

# phase 4 여기서 ,mtx size가 문제가 됨. 모든 행렬을 다 array로 바꾸던지 Matrix로 바꾸던지 해야할 듯.
x_pre = np.array(x_val, dtype=float)
print("x_pre : ", x_pre)
print("del_x : ", del_x)
x_val = x_pre + del_x  # x_val 갱신
print(x_val)

# 4단계: x_pre 갱신
# x_pre = np.array(x_val, dtype=float)  # x_val을 numpy array로 변환
# print("x_pre:", x_pre)
# x_val = x_pre + del_x  # x_val 갱신
# print("x_val 갱신 결과:", x_val)