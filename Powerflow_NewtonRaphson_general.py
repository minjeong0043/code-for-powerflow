import pandas as pd
import numpy as np
import sympy as sp
import time
from Powerflow_NewtonRaphson_general_func import Ybus,Ybus_with_transformer, DefineVariable, Ufunc, Jacobian, PQ_given, NewtonRaphson, Result_values, NewtonRaphson2, NewtonRaphson3

data_name = "ac_case5.xlsx"
# data_name = "ac_case3_fin.xlsx"
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
if transformer.size == 0:
    Y_mag, Y_rad = Ybus(bus_pu, branch)
else:
    Y_mag, Y_rad, Y = Ybus_with_transformer(bus_pu, branch, transformer)
# print(f"Y_mag : {Y_mag}\nY_rad : {Y_rad}")
# print(f"Y : {Y}")

## Define Variable and initial value
V_mag, delta, sym_variables, indices_delta = DefineVariable(bus_pu, generator_pu)
# print(len(sym_variables))
# print(sym_variables)

## PQ func
Ufunc = Ufunc(bus_pu, Y_mag, Y_rad, V_mag, delta)
# print(f"Ufunc : {Ufunc},  size :{Ufunc.shape}")

## Jacobian
J = Jacobian(bus_pu, Ufunc, sym_variables)
print(f"J : {J}")
## Given Data
U_given = PQ_given(bus_pu, generator_pu)
print(f"U_given :{U_given}")


## iteration
tolerance = 1e-6
# tolerance = 1
max_iter = 10
iteration = 0
converged = False

# 맨 초기값인데 엑셀 데이터 고려해서 수정할 필요 있음.
x = {}
for var in sym_variables:
    if "delta" in str(var):
        x[var] = 0
    if "V" in str(var):
        x[var] = 1
start_time = time.time()
print(f"Jacobian_init : {J.subs(x).evalf()}")
while not converged or iteration < max_iter:
    # print(f"x :{x}")
    # print(f"iter = {iteration}")
    Ufunc_cal = Ufunc.subs(x)
    del_U = np.zeros(len(bus_pu))
    # print(f"Ufunc_cal :{Ufunc_cal}  {type(Ufunc_cal)}, del_U :{del_U}")
    # print(f"Ugiven :{U_given},  {type(U_given)}")
    for i in range(len(bus_pu)):
        del_U[i] = U_given[i] - Ufunc_cal[i]

    print(f"abs(del_U) :{abs(del_U)}")
    if np.all(abs(del_U) < tolerance):
    # if np.all(abs(del_x) < tolerance):
        converged = True
        print(f"converged at i = {iteration}, x = {x}, Ufunc_cal = {Ufunc_cal}")
        break
    x, del_x = NewtonRaphson3(x, Ufunc, bus_pu, U_given, J, indices_delta)
    iteration += 1
    print(f"iter = {iteration}")
    print("------------- %s seconds--------------" % (time.time() - start_time))


## Results
V_result, delta_result, P_result, Q_result = Result_values(x, bus_pu, V_mag, delta, Y_mag, Y_rad)

bus = []
delta_result_degree = np.zeros(len(delta_result))
for i in range(len(bus_pu)):
    bus.append(i+1)
    delta_result_degree = delta_result * 180 / np.pi

## Save Data to excel file
P_G = np.zeros(len(bus_pu))
Q_G = np.zeros(len(bus_pu))
P_load = np.zeros(len(bus_pu))
Q_load = np.zeros(len(bus_pu))
for i in range(len(generator_pu)):
    print(generator_pu['Bus'] - 1)
    P_G[generator_pu['Bus'] - 1] = generator_pu['PG (pu)'][i]
    Q_G[generator_pu['Bus'] - 1] = generator_pu['QG (pu)'][i]
for i in range(len(bus_pu)):
    P_load[i] = bus_pu['Pload (pu)'][i]
    Q_load[i] = bus_pu['Qload (pu)'][i]

for i in range(len(bus_pu)):
    if bus_pu['Type'][i] == 'Swing':
        P_G[i] = P_result[i]
        Q_G[i] = Q_result[i]
    elif bus_pu['Type'][i] == 'PV':
        Q_G[i] = Q_result[i] + Q_load[i]

data = {
    'Bus' : bus,
    'Volt': V_result,
    'Angle': delta_result_degree,
    'PG': P_G,
    'QG': Q_G,
    'Pload': P_load,
    'Qload': Q_load
}

df = pd.DataFrame(data)
df.to_excel('powerflow_results.xlsx', index=False)

print("The file has been saved!!!")