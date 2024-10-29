import pandas as pd
import numpy as np
import sympy as sp
from Powerflow_NewtonRaphson_general_func import Ybus,Ybus_with_transformer, DefineVariable, Ufunc, Jacobian, PQ_given, NewtonRaphson, Result_values

data_name = "ac_case3_fin.xlsx"
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
    Y_mag, Y_rad = Ybus_with_transformer(bus_pu, branch, transformer)

## Define Variable and initial value
V_mag, delta, sym_variables, indices_delta = DefineVariable(bus_pu, generator_pu)

## PQ func
Ufunc = Ufunc(bus_pu, Y_mag, Y_rad, V_mag, delta)

## Jacobian
J = Jacobian(bus_pu, Ufunc, sym_variables)

## Given Data
U_given = PQ_given(bus_pu, generator_pu)


## iteration
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
    # print(f"x :{x}")
    Ufunc_cal = Ufunc.subs(x)
    del_U = np.zeros(len(bus_pu))
    # print(f"Ufunc_cal :{Ufunc_cal}, del_U :{del_U}")
    for i in range(len(bus_pu)):
        del_U[i] = U_given[i] - Ufunc_cal[i]
    if np.all(del_U < tolerance):
        converged = True
        print(f"converged at i = {iter},\n x = {x},\n Ufunc_cal = {Ufunc_cal}")
        break
    x, x_degree = NewtonRaphson(x, Ufunc, bus_pu, U_given, J, indices_delta)
    iter += 1
    # print(f"iter = {iter}")

# Results
V_result, delta_result, P_result, Q_result = Result_values(x, bus_pu, V_mag, delta, Y_mag, Y_rad)

# Data to excel
data = {

}