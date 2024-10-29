import pandas as pd
import numpy as np
import sympy as sp
from Powerflow_NewtonRaphson_func import Ybus, DefineVariable, Ufunc, Jacobian, PQ_given

data_name = "ac_case3.xlsx"

bus = pd.read_excel(data_name, sheet_name='bus')
branch = pd.read_excel(data_name, sheet_name='branch')

Y_mag, Y_rad = Ybus(bus, branch)
# print(f"Y_mag : {Y_mag},       Y_rad : {Y_rad}")

V_mag, delta, sym_variables, indices_delta = DefineVariable(bus)
# print(f"V_mag : {V_mag},     delta : {delta},         indices_delta : {indices_delta}") #delta의 단위가 무엇이냐

Ufunc = Ufunc(bus, Y_mag, Y_rad, V_mag, delta)
# print(f"Ufunc : {Ufunc}")

J = Jacobian(bus, Ufunc, sym_variables)
# print(f"J : {J}")

# Given Data
U_given = PQ_given(bus)
# print(f"U_given : {U_given}")

# iteration 1
var0_dict = {}
for var in sym_variables:
    if "delta" in str(var):
        var0_dict[var] = 0
    if "V" in str(var):
        var0_dict[var] = 1

Ufunc0_cal = Ufunc.subs(var0_dict)
del_U = np.zeros(len(bus))
for i in range(len(bus)):
    del_U[i] = U_given[i] - Ufunc0_cal[i]
del_U = sp.Matrix(del_U)
# print(f"del_U : {del_U}")
J0 = J.subs(var0_dict)
# print(f"J0 : {J0}")
del_x = J0.inv() @ del_U
# print(del_x)

# rad to degree
del_x_degree = del_x.copy()
for i in range(len(indices_delta)):
    del_x_degree[i] = del_x[i] * 180 / np.pi
    # print(del_x[i] * 180 / np.pi)
# print(f"del_x_degree : {del_x_degree}")

# Update x
x1 = {}
x0 = var0_dict.values()
x0_values = list(x0)
for i, key in enumerate(var0_dict):
    x1[key] = np.array(x0_values[i]) + del_x[i]
# print(f"x1 : {x1}")

# iteration 2
Ufunc1_cal = Ufunc.subs(x1)
# print(f"Ufunc1_cal = {Ufunc1_cal}")
del_U = np.zeros(len(bus))
for i in range(len(bus)):
    del_U[i] = U_given[i] - Ufunc1_cal[i]
del_U = sp.Matrix(del_U)
# print(f"del_U : {del_U}")
J1 = J.subs(x1)
# print(f"J1 : {J1}")
del_x = J1.inv() @ del_U
# print(f"del_x : {del_x}")
del_x_degree = del_x.copy()
for i in range(len(indices_delta)):
    del_x_degree[i] = del_x[i] * 180/ np.pi
# print(f"del_x_degree : {del_x_degree}")

# Update x
x2 = {}
x2_degree = {}
x1 = x1.values()
x1_values = list(x1)
for i, key in enumerate(var0_dict):
    x2[key] = np.array(x1_values[i]) + del_x[i]
    if 'delta' in str(key):  # key에 'delta'가 포함된 경우
        x2_degree[key] = (np.array(x1_values[i]) + del_x[i]) * 180 / np.pi  # degree 변환
    else:
        x2_degree[key] = np.array(x1_values[i]) + del_x[i]  # degree 변환 없이 그대로 저장

# print(f"x2 : {x2}\nx2_degree : {x2_degree}")


# iteration 3
Ufunc2_cal = Ufunc.subs(x2)
# print(f"Ufunc2_cal = {Ufunc2_cal}")
del_U = np.zeros(len(bus))
for i in range(len(bus)):
    del_U[i] = U_given[i] - Ufunc2_cal[i]
del_U = sp.Matrix(del_U)
# print(f"del_U : {del_U}")
J2 = J.subs(x2)
del_x = J2.inv() @ del_U
# print(f"del_x : {del_x}")
del_x_degree = del_x.copy()
for i in range(len(indices_delta)):
    del_x_degree[i] = del_x[i] *180/ np.pi
# print(f"del_x_degree : {del_x_degree}")

# Update x
x3 = {}
x3_degree = {}
x2 = x2.values()
x2_values = list(x2)
for i, key in enumerate(var0_dict):
    x3[key] = np.array(x2_values[i]) + del_x[i]
    if 'delta' in str(key):
        x3_degree[key] = (np.array(x2_values[i]) + del_x[i]) *180/np.pi
    else:
        x3_degree[key] = np.array(x2_values[i]) + del_x[i]

# print(f"x3 :{x3}\nx3_degree : {x3_degree}")

Ufunc3_cal = Ufunc.subs(x3)
print(f"Ufunc3_cal = {Ufunc3_cal}")

# 정해진 변수 기반으로 슬랙 모선에서의 PQ. PV모선에서의 Q값계산
print(f"x3 :{x3}\nx3_degree : {x3_degree}")