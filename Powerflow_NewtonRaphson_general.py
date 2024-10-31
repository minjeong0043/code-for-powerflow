import pandas as pd
import numpy as np
import sympy as sp
from Powerflow_NewtonRaphson_general_func2 import Ybus, Ybus_with_transformer, DefineVariable, Ufunc, Jacobian, PQ_given, NewtonRaphson

data_name = "ac_case5.xlsx"

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

# Cal Ybus
# ---------------------------------------------------------
# Y             : Admittance matrix.
# Y_mag         : Magnitude of the admittance.
# Y_rad         : Admittance angle in radians.
# Y_degrees     : Admittance angle in degree.
# ---------------------------------------------------------
if transformer.size == 0:
    Y_mag, Y_rad, Y_degrees, Y = Ybus(bus_pu, branch)
else:
    Y_mag, Y_rad, Y_degrees, Y = Ybus_with_transformer(bus_pu, branch, transformer)
# print(f"Y : {Y}") # ok
print('\n')

# Define Variable and initial Value
# ---------------------------------------------------------
# sym_variable  : Array of symbolic variable names [delta, V].
# V_mag         : Magnitudes of voltage (V) or symbolic variables, ordered by bus sequence.
# delta         : Phase angle values or symbolic variables, ordered by bus sequence.
# indices_delta : Index positions of delta variables within sym_variable.
# ---------------------------------------------------------
V_mag, delta, sym_variables,sym_variables_init, indices_delta = DefineVariable(bus_pu, generator_pu)
# print(f"V_mag               : {V_mag}")
# print(f"delta               : {delta}")
# print(f"sym_variable        : {sym_variables}")
# print(f"sym_variable_init   : {sym_variables_init}")
# print(f"len(indices_delta)  : {len(indices_delta)}")
# print('\n')


# Bus     Volt    Angle     Real   Reactive
# 1.0000   1.0600        0   1.3122   0.9759
# 2.0000   1.0000  -2.0658   0.2000  -0.5617
# 3.0000   0.9822  -4.5607  -0.4500  -0.1500
# 4.0000   0.9789  -4.8791  -0.4000  -0.0500
# 5.0000   0.9675  -5.7104  -0.6000  -0.1000

# PG func
# ---------------------------------------------------------
# func          : Represents the P and Q functions for each bus, where the index corresponds to the bus number
#                 (PQ bus: both P and Q functions, PV bus: P function only, Swing bus: None).
# Ufunc         : Rearranges the func into [P, Q] format, with P and Q functions ordered sequentially by bus number
#                 (i.e., P1, P2, P3,... followed by Q1, Q2, Q3,...).
# U_given       : Fixed values for the Ufunc.
# P_total       : Total active power, defined as the difference between generated power (PG) and load power (PL).
# Q_total       : Total reactive power, defined as the difference between generated reactive power (QG) and load reactive power (QL).
# ---------------------------------------------------------
func, Ufunc, Ufunc_array_name= Ufunc(bus_pu, Y_mag, Y_rad, V_mag, delta)
P_total, Q_total, U_given = PQ_given(bus_pu, generator_pu)
# print(f"func                : {func}")
# print(f"Ufunc_array_name    : {Ufunc_array_name}")
# print(f"Ufunc               : {Ufunc}")
# print(f"U_given             : {U_given}")
# print('\n')

# Jacobian
# ---------------------------------------------------------
# J             : The Jacobian matrix represented as a NumPy array.
# J_matrix      : To facilitate substitution of values from x, convert J from a NumPy array to a SymPy matrix (sp.Matrix).
# x_init        : A dictionary containing the initial values(sym_variables_init) for each variable(sym_variables).
# J_init        : Initial Jacobian matrix (J) with x_init values substituted.
# ----------------------------------------------------------
J = Jacobian(bus_pu, Ufunc, sym_variables)
J_matrix = sp.Matrix(J)
x_init = {}
for i in range(len(sym_variables)):
    x_init[sym_variables[i]] = sym_variables_init[i]
J_init = np.array(J_matrix.subs(x_init)).astype(np.float64)
# print(f"x_init              : {x_init}")
# print(f"J {J.shape}        : \n{J}")
# print(f"J_init              : \n{np.array(J_matrix.subs(x_init)).astype(np.float64)}") # 초기 자코비안 값이 받은 값과 약간의 차이가 남..! 확인 필요
# print('\n')

# iteration
tolerance = 1e-6
max_iter = 100
iteration = 0
converged = False

x = x_init
J = J_init
# print(f"x                   : {x}")
U_cal = Ufunc.subs(x)
U_given = sp.Matrix(U_given)
# del_U = U_given - U_cal
# print(f"U_given{type(U_given)}: {U_given}")
# print(f"U_cal{type(U_cal)} : {U_cal}")
# print(f"del_U              : {del_U}")
# if all([abs(element) < 1.38 for element in del_U]):
#     print("FFF")


while not converged and iteration < max_iter:
    U_cal = Ufunc.subs(x)
    del_U = U_given - U_cal
    if all([abs(element) < tolerance for element in del_U]):
        converged = True
        print(f"Converged at i = {iteration}, x = {x}, U_cal = {U_cal}")
        break
    x, del_x = NewtonRaphson(x, J_matrix, del_U, indices_delta)
    iteration += 1
    print(f"i = {iteration}")

# Result
# V_mag 와 delta사용하면 될 듯
V_mag_result = [x.get(var, var) for var in V_mag]
delta_result = [x.get(var, var) for var in delta]
print(f"V_mag_result       : {V_mag_result}")
print(f"delta_result       : {delta_result}")

# =======
'''
# Initialize P and Q arrays
n = len(bus_pu)
P = np.zeros(n)
Q = np.zeros(n)

# Power flow calculations for each bus
for i in range(n):
    for j in range(n):
        P[i] += V_mag_result[i] * V_mag_result[j] * Y_mag[i, j] * np.cos(delta_result[i] - delta_result[j] - Y_rad[i, j])
        Q[i] += V_mag_result[i] * V_mag_result[j] * Y_mag[i, j] * np.sin(delta_result[i] - delta_result[j] - Y_rad[i, j])

# Print results
print("P values: ", P)
print("Q values: ", Q)

print(f"---------------")
print(np.cos(np.pi))
print(np.cos(180))'''
