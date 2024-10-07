import numpy as np
from sympy import symbols, Eq, solve, Matrix
# 뉴튼랩슨 방법 (1.역함수)
x1, x2 = symbols('x1 x2')


'''f = [0, 0]
f[0] = x1 + x2
f[1] = x1*x2'''

f = Matrix([x1+x2, x1*x2])

f_prime = f[0].diff(x1)
func_num = 2

x = [x1, x2]

'''for i in range(func_num):
    print("i", i)
    for j in range(func_num):
        print("j", j)
        J[i][j] = f[i].diff(x[j])'''
J = Matrix(2, 2, lambda i, j: f[i].diff(x[j]))

J_inv = J.inv()

#Jinv 구했고 다음으로 x 계산!

x_init = Matrix([4, 9])
x_pre = x_init
y = Matrix([15, 50])
x = x_pre + J_inv @ (y - f(x_pre))

