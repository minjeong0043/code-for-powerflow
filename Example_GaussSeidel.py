import pandas as pd
import numpy as np
from Example_Powerflow_func import Ybus, Cal_PQ, GaussSeidel_PQ

# Read Data
data_name = "ac_case25_example.xlsx"
bus = pd.read_excel(data_name, sheet_name='bus')
branch = pd.read_excel(data_name, sheet_name='branch')

# Calc Ybus
size_bus = len(bus)
Y = Ybus(branch, size_bus)

# Initial voltage assignments for all buses
V_value = np.ones(size_bus) # |V|
delta = np.zeros(size_bus) # 위상
V = np.ones(size_bus, dtype=complex)

Swing_bus_index = bus[bus['Type'] == 'Swing'].index[0]
PV_buses = bus[bus['Type'] == 'PV'].index[0]
V_value[Swing_bus_index] = bus['V (pu)'][Swing_bus_index]
V_value[PV_buses] = bus['V (pu)'][PV_buses]

V = V_value * np.exp(1j*delta)

# 버스에 따라 각 가우스자이델 반복법 적용!
# index_without_Swing = np.where(bus['Type'] != 'Swing')[0]

# for bus_index in index_without_Swing:
#     # print(bus['Type'][bus_index])
#     if bus['Type'][bus_index] == 'PQ':
#         print("bus_index : ",bus_index)
#         V = GaussSeidel_PQ(bus, bus_index, V, Y)
#         print(f"V : {V}")
#     elif bus['Type'][bus_index] == 'PV':
#         print("PV")

# PV bus

tolerance =1e-6
max_iter = 100
iteration = 0
converged = False

index_PV = np.where(bus['Type'] == "PV")[0][0]
print("Type(index_PV) : ", type(index_PV))
print(index_PV)

# del_p와 del_V_value
P_given = bus['PG (pu)'] - bus['PL (pu)']
V_given = bus['V (pu)'][index_PV]

# P_cal, _ = Cal_PQ(V, Y, index_PV, len(bus))
# V_cal =  # 고려 없이 계산해보고 비교하기
# del_P = abs(P_given[index_PV] - P_cal)
# del_V = abs(V_given - V_cal)
# print(len(V_given))
# print(V_given-1)
Q = 0 # PV모선에서의 변수
V_given = V[index_PV]
while not converged and iteration < max_iter:
    P_cal, _ = Cal_PQ(V, Y, index_PV, len(bus))
    print("V_given : ", V_given)
    del_P = abs(P_given[index_PV] - P_cal)
    # np.power(del_V, 2) = np.power(V_given, 2) - np.power(abs(V[index_PV]), 2) # del_V 는 크기임. 위상 x, V_given은 크기로 주어졌고 ,,,
    del_V = (np.power(V_given, 2) - np.power(abs(V[index_PV], 2)))**0.5
    # del_P = abs(V)
    print("V[index_PV]                          : ",V[index_PV])

    if ((del_P < tolerance) and (del_V < tolerance)):
        converged = True
        print(f"Converged at i = {iteration}, V = {V[index_PV]}")
        break

    b = 0
    for j in range(len(bus)):
        if index_PV != j:
            b += Y[index_PV, j] * V[j]
    V[index_PV] = (1/Y[index_PV, index_PV]) * (((P_given[index_PV] - 1j* Q) / V[index_PV].conjugate()) - b)
    print(del_P)
    iteration += 1

