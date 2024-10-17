import pandas as pd
import numpy as np
from Example_Powerflow_func import Ybus

data_name = "ac_case25_example.xlsx"
bus = pd.read_excel(data_name, sheet_name='bus')
branch = pd.read_excel(data_name, sheet_name='branch')

# Calc Ybus
size_bus = len(bus)
Y = Ybus(branch, size_bus)

# 각 버스의 전압 크기와 위상 초기화
V = np.ones(size_bus) # |V|
delta = np.zeros(size_bus) # 위상

Swing_bus_index = bus[bus['Type'] == 'Swing'].index[0]
pv_buses = bus[bus['Type'] == 'PV'].index

V[]
# PQ bus
# bus 2
P2_cal = 0
Q2_cal = 0
t = 0
for i in range(size_bus):
    V