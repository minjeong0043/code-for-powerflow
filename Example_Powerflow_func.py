import numpy as np
import pandas as pd
import cmath

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

# def GaussSeidel_PQ():
#     while True:
#         del =

def Cal_PQ(V, Y, i, size_bus): # 전압 크기, 위상, 어드미턴스, bus i에서의 계산 PQ 전력
    # V_value, delta = cmath.polar(V)
    V_value = []
    delta = []

    # V 배열의 각 요소에 대해 크기와 각도 계산
    for v in V:
        v_val, v_delta = cmath.polar(v)
        V_value.append(v_val)
        delta.append(v_delta)

    # print("delta[i]", delta[i])

    P_cal = 0
    Q_cal = 0
    # t = 0
    for j in range(size_bus):
        P_cal += V_value[i] * V_value[j] * (Y[i][j].real * np.cos(delta[i] - delta[j]) + Y[i][j].imag * np.sin(delta[i] - delta[j]))
        Q_cal += V_value[i] * V_value[j] * (Y[i][j].real * np.sin(delta[i] - delta[j]) - Y[i][j].imag * np.cos(delta[i] - delta[j]))

    # print("P_cal : ", P_cal)
    # print("Q_cal : ", Q_cal)
    return P_cal, Q_cal

def GaussSeidel_PQ(bus, bus_index, V, Y):
    P_given = bus['PG (pu)'] - bus['PL (pu)']
    Q_given = bus['QG (pu)'] - bus['QL (pu)']

    tolerance = 1e-6
    max_iter = 100
    iteration = 0
    converged = False

    while not converged and iteration < max_iter:
        P_cal, Q_cal = Cal_PQ(V, Y, bus_index, len(bus))
        del_P = abs(P_given[bus_index] - P_cal)
        del_Q = abs(Q_given[bus_index] - Q_cal)
        # print(f"del_P : {del_P},                          del_Q : {del_Q},                 V : {V[bus_index]}")
        if ((del_P < tolerance) and (del_Q < tolerance)):
            converged = True
            print(f"Converged at i = {iteration}, V = {V[bus_index]}")
            break

        b = 0
        for j in range(len(bus)):
            if bus_index != j:
                b += Y[bus_index,j] * V[j]
        V[bus_index] = (1/Y[bus_index, bus_index]) * (((P_given[bus_index] - 1j* Q_given[bus_index]) / V[bus_index].conjugate()) - b)

        iteration += 1
        # print(f"iteration = {iteration}")

    return V

def GaussSeidel_PV(bus, bus_index, V, Y):
    P_given = bus['PG (pu)'] - bus['PL (pu)']
    V_given = bus['V (pu)'][bus_index]
    Q = 0
    V_given = V[bus_index]

    tolerance = 1e-6
    max_iter = 100
    iteration = 0
    converged = False

    while not converged and iteration < max_iter:
        P_cal, Q_cal = Cal_PQ(V, Y, bus_index, len(bus))
        del_P = abs(P_given[bus_index] - P_cal)
        del_V = abs(abs(V_given) - abs(V[bus_index]))
        Q = Q_cal
        # print(f"del_P : {del_P},                          del_V : {del_V},                 V : {V[bus_index]}")
        if ((del_P < tolerance) and (del_V < tolerance)):
            converged = True
            print(f"Converged at i = {iteration}, V = {V[bus_index]}")
            break

        b = 0
        for j in range(len(bus)):
            if bus_index != j:
                b += Y[bus_index,j] * V[j]
        V_new = (1/Y[bus_index, bus_index]) * (((P_given[bus_index] - 1j* Q) / V[bus_index].conjugate()) - b)
        #위 식 다른 거로도 돌려보고 결과 비교해보기
        V[bus_index] = V_given * np.exp(1j*np.angle(V_new))
        iteration += 1
        # print(f"iteration = {iteration}")

    return V
