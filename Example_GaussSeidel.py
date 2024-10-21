import pandas as pd
import numpy as np
from Example_Powerflow_func import Ybus, Cal_PQ, GaussSeidel_PQ, GaussSeidel_PV, Ybus_considering_trans

# Read Data
data_name = "ac_case25_example.xlsx"
bus = pd.read_excel(data_name, sheet_name='bus')
branch = pd.read_excel(data_name, sheet_name='branch')
trans = pd.read_excel(data_name, sheet_name='transformer')


# Calc Ybus
Y = Ybus_considering_trans(bus, branch, trans)
print(Y)

# Initial voltage assignments for all buses
size_bus = len(bus)
V_value = np.ones(size_bus) # |V|
delta = np.zeros(size_bus) # 위상
V = np.ones(size_bus, dtype=complex)

Swing_bus_index = bus[bus['Type'] == 'Swing'].index[0]
PV_buses = bus[bus['Type'] == 'PV'].index[0]
V_value[Swing_bus_index] = bus['V (pu)'][Swing_bus_index]
V_value[PV_buses] = bus['V (pu)'][PV_buses]
V = V_value * np.exp(1j*delta)


# 버스에 따라 각 가우스 자이델 반복법 적용!
index_without_Swing = np.where(bus['Type'] != 'Swing')[0]
for bus_index in index_without_Swing:
    # print(bus['Type'][bus_index])
    if bus['Type'][bus_index] == 'PQ':
        print(f"bus_index : {bus_index},          {bus['Type'][bus_index]}")
        V = GaussSeidel_PQ(bus, bus_index, V, Y)
        print(f"V : {V}")
    elif bus['Type'][bus_index] == 'PV':
        print(f"bus_index : {bus_index},          {bus['Type'][bus_index]}")
        V = GaussSeidel_PV(bus, bus_index, V, Y)
        print(f"V : {V}")

