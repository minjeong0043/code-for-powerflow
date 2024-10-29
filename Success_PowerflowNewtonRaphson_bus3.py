import pandas as pd
import numpy as np
import sympy as sp
from Success_PowerflowNewtonRaphson_bus3_func import Ybus, DefineVariable, Ufunc, Jacobian, PQ_given, NewtonRaphson

data_name = "ac_case3.xlsx"

bus = pd.read_excel(data_name, sheet_name='bus')
branch = pd.read_excel(data_name, sheet_name='branch')

# Cal Y bus
Y_mag, Y_rad = Ybus(bus, branch)

# Define Variable and initial value
V_mag, delta, sym_variables, indices_delta = DefineVariable(bus)

# PQ func
Ufunc = Ufunc(bus, Y_mag, Y_rad, V_mag, delta)

# Jacobian
J = Jacobian(bus, Ufunc, sym_variables)

# Given Data
U_given = PQ_given(bus)

# iteration
tolerance = 1e-6
max_iter = 100
iter = 0
converged = False

x = {}
for var in sym_variables:
    if "delta" in str(var):
        x[var] = 0
    if "V" in str(var):
        x[var] = 1
while not converged and iter < max_iter:
    Ufunc_cal = Ufunc.subs(x)
    del_U = np.zeros(len(bus))
    for i in range(len(bus)):
        del_U[i] = U_given[i] - Ufunc_cal[i]
    if np.all(del_U < tolerance):
        converged = True
        print(f"conveged at i = {iter}, x = {x}, Ufunc_cal = {Ufunc_cal}")
        break
    x, x_degree = NewtonRaphson(x, Ufunc, bus, U_given, J, indices_delta)
    iter += 1
    # print(f"iter = {iter}")

# results
V_result = np.empty(len(bus), dtype=float)
delta_result = np.empty(len(bus), dtype=float)
for i in range(len(bus)):
    if isinstance(V_mag[i], np.float64):
        V_result[i] = V_mag[i]
    else:
        V_result[i] = x.get(V_mag[i], np.nan)
    if isinstance(delta[i], np.float64):
        delta_result[i] = delta[i]
    else:
        delta_result[i] = x.get(delta[i], np.nan)
print(f"V_result : {V_result},         delta_result : {delta_result}")

# Cal P, Q each bus
P_result = np.zeros(len(bus), dtype=float)
Q_result = np.zeros(len(bus), dtype=float)

# P, Q 계산
for i in range(len(bus)):
    if bus['Type'][i] == 'PQ':
        for k in range(len(bus)):
            angle_diff = float(delta_result[i] - delta_result[k] - Y_rad[i, k])
            P_result[i] += Y_mag[i, k] * V_result[k] * np.cos(angle_diff)
            Q_result[i] += Y_mag[i, k] * V_result[k] * np.sin(angle_diff)
        P_result[i] = V_result[i] * P_result[i]
        Q_result[i] = V_result[i] * Q_result[i]

    elif bus['Type'][i] == 'PV':
        for k in range(len(bus)):
            angle_diff = float(delta_result[i] - delta_result[k] - Y_rad[i, k])
            P_result[i] += Y_mag[i, k] * V_result[k] * np.cos(angle_diff)
            Q_result[i] += Y_mag[i, k] * V_result[k] * np.sin(angle_diff)
        P_result[i] = V_result[i] * P_result[i]
        Q_result[i] = V_result[i] * Q_result[i]

    elif bus['Type'][i] == 'Swing':
        for k in range(len(bus)):
            angle_diff = float(delta_result[i] - delta_result[k] - Y_rad[i, k])
            P_result[i] += Y_mag[i, k] * V_result[k] * np.cos(angle_diff)
            Q_result[i] += Y_mag[i, k] * V_result[k] * np.sin(angle_diff)
        P_result[i] = V_result[i] * P_result[i]
        Q_result[i] = V_result[i] * Q_result[i]

# 결과 출력
print(f"P_result : {P_result}")
print(f"Q_result : {Q_result}")