import numpy as np
import pandas as pd

def Ybus(branch, size_Mtx):
    Y = np.zeros([size_Mtx, size_Mtx], dtype= complex)

    for i in range(size_Mtx):
        for j in range(size_Mtx):
            if i == j:
                index = np.where((branch['From'] ==i + 1) | (branch['To'] == i + 1))[0]
                if index.size== 0:
                    Y[i,j] = 0
                else:
                    for k in index:
                        Y[i,j] += 1/(branch['R (pu)'][k] + 1j*branch['X (pu)'][k]) + 1j*branch['B (pu)'][k]/2
            else:
                index = np.where(((branch['From'] == i + 1) & (branch['To'] == j + 1)) | ((branch['From'] == j + 1) & (branch['To'] == i + 1)))[0]
                if index.size == 0:
                    Y[i,j] = 0
                else:
                    for k in index:
                        Y[i,j] -= 1/(branch['R (pu)'][k] + 1j*branch['X (pu)'][k])

    return Y