import pandas as pd
import numpy as np
import sympy as sp
from PowerflowNewtonRaphson_bus25_func import Ybus, DefineVariable, Ufunc, Jacobian, PQ_given, NewtonRaphson

data_name = "ac_case25.xlsx"
# Read Data
bus = pd.read_excel(data_name, sheet_name='bus')
branch = pd.read_excel(data_name, sheet_name='branch')
transformer = pd.read_excel(data_name, sheet_name='transformer')
generator = pd.read_excel(data_name, sheet_name='generator')
Sbase = pd.read_excel(data_name, sheet_name='param', header = None)[1][0]

# Change Unit from MVAR, MW to pu
bus_pu = bus.copy()
bus_pu['Pload (MW)'] = bus['Pload (MW)'] / Sbase
bus_pu['Qload (MVAR)'] = bus['Qload (MVAR)'] / Sbase
bus_pu = bus_pu.rename(columns = {'Pload (MW)' : 'Pload (pu)', 'Qload (MVAR)':'Qload (pu)'})

generator_pu = generator.copy()
generator_pu['PG (MW)'] = generator['PG (MW)'] / generator['MBASE (MW)']
generator_pu['QG (MVAR)'] = generator['QG (MVAR)'] / generator['MBASE (MW)']
generator_pu = generator_pu.rename(columns = {'PG (MW)' : 'PG (pu)', 'QG (MVAR)' : 'QG (pu)'})

## Cal Ybus
# Ybus considering Transformer
Y_mag, Y_rad = Ybus(bus_pu, branch, transformer)
# print(f"len(Y_mag) : {len(Y_mag)}")

## Define Variable and initial value
V_mag, delta, sym_variables, indices_delta = DefineVariable(bus_pu, generator_pu)
# print(V_mag)
# print(delta)

## PQ func
Ufunc = Ufunc(bus_pu, Y_mag, Y_rad, V_mag, delta)

## Jacobian
J = Jacobian(bus_pu, Ufunc, sym_variables)
# print(f"len(J) : {len(J)}\nJ : {J}")

## Given Data
U_given = PQ_given(bus_pu, generator_pu)
# print(f"U_given : {U_given}")

## iteration
print("Enter iteration")
tolerance = 1e-6
max_iter = 100
iter = 0
converged = False

# 맨 초기값인데 엑셀 데이터 고려해서 수정할 필요 있음.
x = {}
for var in sym_variables:
    if "delta" in str(var):
        x[var] = 0
    if "V" in str(var):
        x[var] = 1
# print(x)
while not converged and iter < max_iter:
    Ufunc_cal = Ufunc.subs(x)
    del_U = np.zeros(len(bus_pu))
    for i in range(len(bus_pu)):
        del_U[i] = U_given[i] - Ufunc_cal[i]
    if np.all(del_U < tolerance):
        converged = True
        print(f"converged at i = {iter},\n x = {x},\n Ufunc_cal = {Ufunc_cal}")
        break
    x, x_degree = NewtonRaphson(x, Ufunc, bus_pu, U_given, J, indices_delta)
    iter += 1
    print(f"iter = {iter}")

# Results
V_result = np.empty(len(bus_pu), dtype=float)
delta_result = np.empty(len(bus_pu), dtype=float)
for i in range(len(bus_pu)):
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
P_result = np.zeros(len(bus_pu), dtype=float)
Q_result = np.zeros(len(bus_pu), dtype=float)

# P, Q 계산
for i in range(len(bus_pu)):
    if bus_pu['Type'][i] == 'PQ':
        for k in range(len(bus_pu)):
            angle_diff = float(delta_result[i] - delta_result[k] - Y_rad[i, k])
            P_result[i] += Y_mag[i, k] * V_result[k] * np.cos(angle_diff)
            Q_result[i] += Y_mag[i, k] * V_result[k] * np.sin(angle_diff)
        P_result[i] = V_result[i] * P_result[i]
        Q_result[i] = V_result[i] * Q_result[i]

    elif bus_pu['Type'][i] == 'PV':
        for k in range(len(bus_pu)):
            angle_diff = float(delta_result[i] - delta_result[k] - Y_rad[i, k])
            P_result[i] += Y_mag[i, k] * V_result[k] * np.cos(angle_diff)
            Q_result[i] += Y_mag[i, k] * V_result[k] * np.sin(angle_diff)
        P_result[i] = V_result[i] * P_result[i]
        Q_result[i] = V_result[i] * Q_result[i]

    elif bus_pu['Type'][i] == 'Swing':
        for k in range(len(bus_pu)):
            angle_diff = float(delta_result[i] - delta_result[k] - Y_rad[i, k])
            P_result[i] += Y_mag[i, k] * V_result[k] * np.cos(angle_diff)
            Q_result[i] += Y_mag[i, k] * V_result[k] * np.sin(angle_diff)
        P_result[i] = V_result[i] * P_result[i]
        Q_result[i] = V_result[i] * Q_result[i]

# 결과 출력
print(f"P_result : {P_result}")
print(f"Q_result : {Q_result}")