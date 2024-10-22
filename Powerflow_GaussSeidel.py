import pandas as pd
import numpy as np
from Powerflow_func import Ybus_considering_trans, GaussSeidel_PQ, GaussSeidel_PV, ReadData_Unit_pu, Ybus

# read input file and define parameters
data_path = "C:\\Users\\leejunwoo\\PycharmProjects\\pythonProject\\ac_case25.xlsx"

bus_pu, branch, transformer, generator_pu, Sbase = ReadData_Unit_pu(data_path)

# Y = Ybus_considering_trans(bus_pu, branch, transformer)
Y = Ybus(branch, len(bus_pu))

# initial voltage assignments for all buses
V_value = bus_pu['Vm (pu)'].copy()
V_angle = bus_pu['Va (degree)'].copy()

# PV버스의 경우엔 SETPOINT로 전압을 바꿔줘야함.bus k의 setpoint
for i in range(len(generator_pu)):
    V_value.loc[generator_pu['Bus'][i] -1] = generator_pu['Voltage setpoint (pu)'][i]
V = V_value * np.exp(1j * V_angle)
# V = np.array(V) # pandas series --> np.array


# Gauss Seidel
index_without_Swing = np.where(bus_pu['Type'] != 'Swing')[0]
# print(index_without_Swing)
for bus_index in index_without_Swing:
    if bus_pu['Type'][bus_index] == 'PQ':
        V = GaussSeidel_PQ(data_path, bus_index, V, Y)
        print(f"bus_index : {bus_index+1} ({bus_pu['Type'][bus_index]}),     V_mag : {abs(V[bus_index])}     V_angle : {np.angle(V[bus_index])}")
    elif bus_pu['Type'][bus_index] == 'PV':
        V = GaussSeidel_PV(data_path, bus_index, V, Y)
        print(f"bus_index : {bus_index+1} ({bus_pu['Type'][bus_index]}),     V_mag : {abs(V[bus_index])}     V_angle : {np.angle(V[bus_index])}")


