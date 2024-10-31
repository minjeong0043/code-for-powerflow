import pandas as pd
import numpy as np
import sympy as sp
from scipy.linalg import lu_factor, lu_solve

def Ybus(bus, branch):
    Y = np.zeros([len(bus), len(bus)], dtype=complex)
    Y_mag = np.zeros([len(bus), len(bus)], dtype=int)
    Y_degrees = np.zeros([len(bus), len(bus)], dtype=int)

    for i in range(len(bus)):
        for j in range(len(bus)):
            if i == j:
                index = np.where((branch['From'] == i + 1) | (branch['To'] == i + 1))[0]
                for k in index:
                    Y[i,j] += 1/(branch['R (pu)'][k] + 1j * branch['X (pu)'][k]) + 1j* branch['B (pu)'][k]/2
                Y_mag[i, j] = np.abs(Y[i, j])
                Y_degrees[i, j] = np.degrees(np.angle(Y[i, j]))
            elif i != j:
                index = np.where(((branch['From'] == i + 1) & (branch['To'] == j + 1)) | ((branch['From'] == j + 1) & (branch['To'] == i + 1)))[0][0]
                if index.size == 0:
                    Y[i, j] = 0
                    Y_mag[i, j] = np.abs(Y[i, j])
                    Y_degrees[i, j] = np.degrees(np.angle(Y[i, j]))
                else:
                    Y[i, j] = - 1/(branch['R (pu)'][index] + 1j * branch['X (pu)'][index])
                    Y_mag[i, j] = np.abs(Y[i, j])
                    Y_degrees[i, j] = np.degrees(np.angle(Y[i, j]))

    Y_rad = np.radians(Y_degrees)
    return Y_mag, Y_rad,Y_degrees, Y

def Ybus_with_transformer(bus, branch, transformer):
    trans = []
    index_trans = []
    for i in range(len(transformer)):
        for j in range(len(branch)):
            if (((branch['From'][j] == transformer['From'][i]) & (branch['To'][j] == transformer['To'][i])) | (
                    (branch['From'][j] == transformer['To'][i]) & (branch['To'][j] == transformer['From'][i]))):
                # From, To, tap, index in branch sheet
                trans.append([transformer['From'][i], transformer['To'][i], transformer['Tap'][i], j])
                index_trans.append(j)

    Y = np.zeros([len(bus), len(bus)], dtype=complex)
    Y_mag = np.zeros([len(bus), len(bus)], dtype=int)
    Y_degrees = np.zeros([len(bus), len(bus)], dtype=int)
    for i in range(len(bus)):
        for j in range(len(bus)):
            if i == j:
                index = np.where((branch['From'] == i + 1) | (branch['To'] == i + 1))[0]
                if len(set(index) & set(index_trans)) == 0:
                    for k in index:
                        Y[i, j] += 1 / (branch['R (pu)'][k] + 1j * branch['X (pu)'][k]) + 1j * branch['B (pu)'][k] / 2
                elif len(set(index) & set(index_trans)) != 0:
                   for k in (set(index) - set(index_trans)):
                        Y[i, j] += 1 / (branch['R (pu)'][k] + 1j * branch['X (pu)'][k]) + 1j * branch['B (pu)'][k] / 2
                   for k in (set(index) & set(index_trans)):
                        t = [item for item in trans if item[3] == k][0][2]
                        Y[i, j] += (1 / (branch['R (pu)'][k] + 1j * branch['X (pu)'][k])) / np.power(t, 2)
            elif i != j:
                index = np.where(((branch['From'] == i + 1) & (branch['To'] == j + 1)) | (
                            (branch['From'] == j + 1) & (branch['To'] == i + 1)))[0]
                if len(set(index) & set(index_trans)) == 0:
                    if index.size == 0:
                        Y[i, j] = 0
                    else:
                        Y[i, j] = -1 / (branch['R (pu)'][index[0]]) + 1j * branch['X (pu)'][index[0]]
                elif len(set(index) & set(index_trans)) != 0:
                    for k in (set(index) - set(index_trans)):
                        Y[i, j] = -1 / (branch['R (pu)'][k]) + 1j * branch['X (pu)'][k]
                    for k in (set(index) & set(index_trans)):
                        t = [item for item in trans if item[3] == k][0][2]
                        Y[i, j] = (-1 / (branch['R (pu)'][k] + 1j * branch['X (pu)'][k])) / t
    for i in range(len(bus)):
        for j in range(len(bus)):
            Y_mag[i, j] = np.abs(Y[i, j])
            Y_degrees[i, j] = np.degrees(np.angle(Y[i, j]))

    Y_rad = np.radians(Y_degrees)

    return Y_mag, Y_rad, Y_degrees, Y

def DefineVariable(bus_pu, generator_pu):
    V_mag_init = bus_pu['Vm (pu)'].copy()
    delta_init = bus_pu['Va (degree)'].copy()
    for i in range(len(generator_pu)):
        V_mag_init.loc[generator_pu['Bus'][i] - 1] = generator_pu['Voltage setpoint (pu)'][i]

    variables = []
    indices_delta = []
    sym_variables_init = []
    for i in range(len(bus_pu)):
        if bus_pu['Type'][i] == 'PQ' or bus_pu['Type'][i] == 'PV':
            variables.append(f"delta{i}")
            sym_variables_init.append(bus_pu['Va (degree)'][i])
            indices_delta.append(i)
    for i in range(len(bus_pu)):
        if bus_pu['Type'][i] == 'PQ':
            variables.append(f"V{i}")
            sym_variables_init.append(bus_pu['Vm (pu)'][i])
    sym_variables = sp.symbols(variables)
    delta = np.empty(len(bus_pu), dtype=object)
    V_mag = np.empty(len(bus_pu), dtype=object)

    for i in range(len(bus_pu)):
        delta_name = f"delta{i}"
        V_name = f"V{i}"
        if delta_name in variables:
            delta[i] = sp.symbols(delta_name)
        else:
            delta[i] = delta_init[i]

        if V_name in variables:
            V_mag[i] = sp.symbols(V_name)
        else:
            V_mag[i] = V_mag_init[i]

    return V_mag, delta, sym_variables,sym_variables_init, indices_delta

def Ufunc(bus, Y_mag, Y_rad, V_mag, delta): # PQ func
    func = np.zeros(len(bus), dtype=object)
    for i in range(len(bus)):
        if bus['Type'][i] == 'PQ':
            P = np.zeros(1, dtype=object)
            Q = np.zeros(1, dtype=object)
            for k in range(len(bus)):
                P += Y_mag[i, k] * V_mag[k] * sp.cos(delta[i] - delta[k] - Y_rad[i, k])
                Q += Y_mag[i, k] * V_mag[k] * sp.sin(delta[i] - delta[k] - Y_rad[i, k])
            P = V_mag[i] * P
            Q = V_mag[i] * Q
            func[i] = (P[0], Q[0])
        elif bus['Type'][i] == 'PV':
            P = np.zeros(1, dtype=object)
            for k in range(len(bus)):
                P += Y_mag[i, k] * V_mag[k] * sp.cos(delta[i] - delta[k] - Y_rad[i, k])
            P = V_mag[i] * P
            func[i] = (P[0],)
        elif bus['Type'][i] == 'Swing':
            func[i] = ()

    func_array = []
    func_array_name = []
    for i in range(len(func)):
        if len(func[i]) != 0:
            if len(func[i]) == 1:
                func_array.append(func[i][0])
                func_array_name.append(f'P{i}')
            elif len(func[i]) == 2:
                func_array.append(func[i][0])
                func_array_name.append(f'P{i}')
    for i in range(len(func)):
        if len(func[i]) == 2:
            func_array.append(func[i][1])
            func_array_name.append(f'Q{i}')
    func_array = sp.Matrix(func_array)

    return func, func_array, func_array_name

def PQ_given(bus_pu, generator_pu):
    PG = np.zeros(len(bus_pu))
    QG = np.zeros(len(bus_pu)) # 사실상 PV모선에서 Q는 변수이기때문에 필요 없을 거라고 생각.. check 필요
    for i in range(len(generator_pu)):
        PG[generator_pu['Bus'][i] -1] = generator_pu['PG (pu)'][i]
        QG[generator_pu['Bus'][i] - 1] = generator_pu['QG (pu)'][i]

    Pload = bus_pu['Pload (pu)'].copy()
    Qload = bus_pu['Qload (pu)'].copy()

    P_total = PG - Pload
    Q_total = QG - Qload

    U_given = []
    for i in range(len(bus_pu)):
        if bus_pu['Type'][i] == 'PQ':
            U_given.append(P_total[i])
        elif bus_pu['Type'][i] == 'PV':
            U_given.append(P_total[i])
    for i in range(len(bus_pu)):
        if bus_pu['Type'][i] == 'PQ':
            U_given.append(Q_total[i])

    return P_total, Q_total, U_given

def Jacobian(bus, func_array, sym_variables):
    J = np.zeros([len(func_array), len(func_array)], dtype=object)
    for i in range(len(func_array)):
        for j in range(len(sym_variables)):
            J[i, j] = func_array[i].diff(sym_variables[j])

    return J

def NewtonRaphson(x, J, del_U, indices_delta):
    J_val = J.subs(x)
    lu, piv = lu_factor(J_val)  # LU 분해
    del_x = lu_solve((lu, piv), del_U)  # LU 분해 결과로 연립 방정식 풀이

    a = enumerate(x)
    x_new = {}
    x_values = list(x.values())
    for i, key in a:
        x_new[key] = np.array(x_values[i]) + del_x[i][0]

    return x_new, del_x

def NewtonRaphson2(x, J, del_U, indices_delta):
    J_val = J.subs(x)
    # print(f"J_val : {J_val}")
    # print(f"del_U : {del_U}")
    # LU 분해 적용 (scipy 라이브러리 사용)
    lu, piv = lu_factor(J_val)  # LU 분해
    del_x = lu_solve((lu, piv), del_U)  # LU 분해 결과로 연립 방정식 풀이
    # print(f"del_x :{del_x}")
    # 새로운 값 계산
    a = enumerate(x)
    x_new = {}
    x_values = list(x.values())

    # 각 변수에 대해 업데이트된 값을 계산
    for i, key in a:
        x_new[key] = np.array(x_values[i]) + del_x[i][0]

    ####
    v_key = [key for key in x if 'V' in str(key)]
    for i in v_key:
        if (x_new[i] < 0.95):
            x_new[i] = 0.95
        elif (x_new[i] > 1.05):
            x_new[i] = 1.05

    return x_new, del_x
