import pandas as pd
import numpy as np
from sympy.printing.precedence import precedence_Integer

# read input file and define parameters
data_path = "ac_case25.xlsx"

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
From = branch['From']
To = branch['To']
R = branch['R (pu)']
X = branch['X (pu)']
B = branch['B (pu)']
size_Mtx = len(bus)

Y = np.zeros([size_Mtx, size_Mtx], dtype = complex)

for i in range(size_Mtx):
    for j in range(size_Mtx):
        if i == j:
            index = np.where(From == i)[0]
            if index.size == 0:
                Y[i,j] = 0
            else:
                for k in index:
                    Y[i,j] += 1 / (R[k] + 1j*X[k]) + 1j*B[k]/2
        else:
            index = np.where(((From == i) & (To == j)) | ((From == j) & (To == i)))[0]
            if index.size == 0:
                Y[i,j] = 0
            else:
                # Y[i,j] = -1/(R[index] + 1j*X[index])
                for k in index:
                    Y[i,j] -= 1/(R[k] + 1j*X[k])


# Gauss Seidel method according to bus Type
Type = bus['Type']
# print(Type)
for i in range(len(Type)):
    if Type[i] == 'PV':
        print("PV bus")
        GaussSeidel_PV()
    elif Type[i] == 'PQ':
        print("PQ bus")
        GaussSeidel_PQ()
    elif Type[i] == 'Swing':
        print("Swing bus")
        GaussSeidel_Swing()
