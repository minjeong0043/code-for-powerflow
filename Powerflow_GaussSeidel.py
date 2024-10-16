import pandas as pd
import numpy as np
# from sympy.printing.precedence import precedence_Integer
from Powerflow_func import Ybus, GaussSeidel_PQ

# read input file and define parameters
data_path = "../../Downloads/ac_case25.xlsx"

bus = pd.read_excel(data_path, sheet_name='bus')
branch = pd.read_excel(data_path, sheet_name='branch')
transformer = pd.read_excel(data_path, sheet_name='transformer')
generator = pd.read_excel(data_path, sheet_name='generator')
Sbase = pd.read_excel(data_path, sheet_name='param', header=None).iloc[0,1]

Pload = bus['Pload (MW)']
Qload = bus['Qload (MVAR)']

# Transform from MW, MVAR to pu
Pload_pu = Pload/Sbase
Qload_pu = Qload/Sbase

# Calculate Ybus
size_Mtx = len(bus)
Y = Ybus(branch, size_Mtx)

# Gauss Seidel method according to bus Type
Type = bus['Type']
# for i in range(len(Type)):
#     if Type[i] == 'PV':
#         print("PV bus")
#         GaussSeidel_PV()
#     elif Type[i] == 'PQ':
#         print("PQ bus")
#         GaussSeidel_PQ()
#     elif Type[i] == 'Swing':
#         print("Swing bus")
#         GaussSeidel_Swing()

# TEST GaussSeidel_PQ()
V = np.ones(size_Mtx, dtype=complex)
for i in range(len(Type)):
    if Type[i] == 'PQ':
        print("i = ", i)
        V = GaussSeidel_PQ(Pload_pu, Qload_pu, Y, i, V)
        print("V : ", V)





