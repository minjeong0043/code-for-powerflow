import numpy as np
import pandas as pd

bus = pd.read_excel("ac_case25.xlsx", sheet_name= 'bus')
branch = pd.read_excel("ac_case25.xlsx", sheet_name= 'branch')
#print(bus)

size_Mtx = len(bus)
# print(size_Mtx)

## Cal Y bus
From = branch['From']
To = branch['To']
R = branch['R (pu)']
X = branch['X (pu)']
B = branch['B (pu)']

Y = np.zeros([size_Mtx, size_Mtx], dtype= complex)

for i in range(size_Mtx):
    for j in range(size_Mtx):
        if i == j:
            index = np.where(From == i)[0]
            if index.size == 0:
                Y[i,j] = 0
            else:
                for k in index:
                    Y[i,j] += 1/ (R[k] + 1j*X[k]) + 1j*B[k]/2
        else:
            index = np.where(((From == i) & (To == j)) |(From == j) & (To == i))[0]
            if index.size == 0:
                Y[i, j] = 0
            else:
                # Y[i,j] = - 1/(R[index] + 1j*X[index])
                for k in index:
                    Y[i, j] -= 1/(R[k] + 1j*X[k])

#print(Y)

Type = bus['Type']

print(Type[0])
print(type(Type[0]))
# Gauss Seidel
# load bus(PQ)
V = np.ones(size_Mtx, dtype=complex)

for i in range(len(Type)): # i is bus num
    if Type[i] == 'PV': # 전압제어모선
        print("PV bus")
    elif Type[i] == 'PQ': # 부하 모선
        print("PQ bus")
        V[i+1] = (1 / Y[1][1]) * (((P[1] - 1j * Q[1]) / V[0].conjugate()) - (Y[1][3] * V[3] + Y[1][4] * V[4]))
        V[1] = (1 / Y[1][1]) * (((P[1] - 1j * Q[1]) / V[1].conjugate()) - (Y[1][3] * V[3] + Y[1][4] * V[4]))
    elif Type[i] == 'Swing': # swing모선
     print("Swing bus")
    else:
        break

# 0 번 버스
if Type[0] == 'PV':
    print("PV bus")
elif Type[0] == 'PQ':
    print("PQ bus")
elif Type[0] == 'Swing':
    print("Swing bus")
