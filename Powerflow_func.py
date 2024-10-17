import numpy as np
def TransformUnit(bus, Sbase):
    bus['Pload (MW)'] = bus['Pload (MW)'] / Sbase
    bus['Qload (MVAR)'] = bus['Qload (MVAR)'] / Sbase
    bus.rname(colimns = {'Pload (MW)' : 'Pload (pu)', 'Qload (MVAR)' : 'Qlaod (pu)'})

    return bus

def GaussSeidel_PQ(P, Q, Y, i, V): #load bus
    size_Mtx = len(P)
    # V = np.ones(size_Mtx, dtype = complex)
    b = np.zeros(1, dtype=complex)
    for s in range(size_Mtx):
        if s!= i:
            # print("s",s)
            b += Y[i][s] *V[s]
            # print("Y : ", Y[i][s])
            # print('V : ', V[s])
            # print(b)
    # iteration 1
    V[i] = (1 / Y[i][i]) * (((P[i] - 1j*Q[i]) / V[i].conjugate()) - b)
    # iteration 2
    V[i] = (1 / Y[i][i]) * (((P[i] - 1j*Q[i]) / V[i].conjugate()) - b)

    return V


# def GaussSeidel_PV():
#
# def GaussSeidel_Swing():

def Ybus(branch, size_Mtx):
    From = branch['From'] -1
    To = branch['To'] -1
    R = branch['R (pu)']
    X = branch['X (pu)']
    B = branch['B (pu)']

    Y = np.zeros([size_Mtx, size_Mtx], dtype= complex)

    for i in range(size_Mtx):
        for j in range(size_Mtx):
            if i == j:
                index = np.where((From ==i) | (To == i))[0]
                if index.size== 0:
                    Y[i,j] = 0
                else:
                    for k in index:
                        Y[i,j] += 1/(R[k] + 1j*X[k]) + 1j*B[k]/2
            else:
                index = np.where(((From == i) & (To == j)) | ((From == j) & (To == i)))[0]
                if index.size == 0:
                    Y[i,j] = 0
                else:
                    for k in index:
                        Y[i,j] -= 1/(R[k] + 1j*X[k])

    return Y