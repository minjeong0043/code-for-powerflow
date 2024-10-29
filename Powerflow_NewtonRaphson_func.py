import pandas as pd
import numpy as np
import sympy as sp

def Ybus(bus, branch):
    Y = np.zeros([len(bus), len(bus)], dtype=complex)
    Y_mag = np.zeros([len(bus), len(bus)], dtype=int)
    Y_degrees = np.zeros([len(bus), len(bus)], dtype=int)

    for i in range(len(bus)):
        for j in range(len(bus)):
            if i == j:
                index = np.where((branch['From'] == i + 1) | (branch['To'] == i + 1))[0]
                # print(f"index : {index}")
                for k in index:
                    # print(f"k : {k}")
                    Y[i, j] += branch['G'][k] + 1j * branch['B'][k]
                Y_mag[i, j] = np.abs(Y[i, j])
                Y_degrees[i, j] = np.degrees(np.angle(Y[i, j]))
            elif i != j:
                index = np.where(((branch['From'] == i + 1) & (branch['To'] == j + 1)) | (
                            (branch['From'] == j + 1) & (branch['To'] == i + 1)))[0][0]
                # print(f"i, j { i+1, j+1}         index : {index}")
                if index.size == 0:
                    Y[i, j] = 0
                    Y_mag[i, j] = np.abs(Y[i, j])
                    Y_degrees[i, j] = np.degrees(np.angle(Y[i, j]))
                else:
                    Y[i, j] = - (branch['G'][index] + 1j * branch['B'][index])
                    Y_mag[i, j] = np.abs(Y[i, j])
                    Y_degrees[i, j] = np.degrees(np.angle(Y[i, j]))
    # print(f"Y_mag : {Y_mag}         Y_degrees: {sp.rad(Y_degrees)}")
    Y_rad = sp.rad(Y_degrees)
    return Y_mag, Y_rad

def DefineVariable(bus):
    V_mag_init = np.ones(len(bus))
    delta_init = np.zeros(len(bus))
    for i in range(len(bus)):
        if bus['Type'][i] == 'Swing':
            V_mag_init[i] = bus['V'][i]
            delta_init[i] = bus['delta'][i]
        elif bus['Type'][i] == 'PV':
            V_mag_init[i] = bus['V'][i]
    # print(f"V: {V_init},       delta : {delta_init}")

    variables = []
    indices_delta = []
    for i in range(len(bus)):  # Add delta variables for PQ and PV buses
        if bus['Type'][i] == 'PQ' or bus['Type'][i] == 'PV':
            variables.append(f"delta[{i}]")
            indices_delta.append(i)

    for i in range(len(bus)):  # Add voltage variables for PQ buses
        if bus['Type'][i] == 'PQ':
            variables.append(f"V[{i}]")
    sym_variables = sp.symbols(variables)
    # print(f"sympy variables: {sym_variables}")

    # 변수와 초기값 동시 가지는 배열
    delta = np.zeros(len(bus), dtype=object)
    V_mag = np.ones(len(bus), dtype=object)

    for i in range(len(bus)):
        delta_name = f"delta[{i}]"
        V_name = f"V[{i}]"
        if delta_name in variables:
            delta[i] = sp.symbols(delta_name)
            # print(f"delta[i] : {delta[i]}")
        else:
            delta[i] = delta_init[i]

        if V_name in variables:
            V_mag[i] = sp.symbols(V_name)
            # print(f"V[i] : {V[i]}")
        else:
            V_mag[i] = V_mag_init[i]
    # print(type(delta[1]))
    # print(f"delta : {delta},          V : {V}")
    return V_mag, delta, sym_variables, indices_delta

def Ufunc(bus, Y_mag, Y_rad, V_mag, delta): # PQ func
    func = np.zeros(len(bus), dtype=object)
    for i in range(len(bus)):
        # print(f"i : {i}")
        if bus['Type'][i] == 'PQ':
            P = np.zeros(1, dtype=object)
            Q = np.zeros(1, dtype=object)
            for k in range(len(bus)):
                P += Y_mag[i, k] * V_mag[k] * sp.cos(delta[i] - delta[k] - Y_rad[i, k])
                Q += Y_mag[i, k] * V_mag[k] * sp.sin(delta[i] - delta[k] - Y_rad[i, k])
            P = V_mag[i] * P
            Q = V_mag[i] * Q
            # print(f"final_P :{P[0]}, \nfinal_Q :{Q}")
            func[i] = (P[0], Q[0])
        elif bus['Type'][i] == 'PV':
            P = np.zeros(1, dtype=object)
            for k in range(len(bus)):
                P += Y_mag[i, k] * V_mag[k] * sp.cos(delta[i] - delta[k] - Y_rad[i, k])
            P = V_mag[i] * P
            # print(f"final_P :{P[0]}")
            func[i] = (P[0],)
        elif bus['Type'][i] == 'Swing':
            func[i] = ()
    # print(f"len(func[0]) : {len(func[2])}")

    func_array = []
    for i in range(len(func)):
        if len(func[i]) != 0:
            if len(func[i]) == 1:
                func_array.append(func[i][0])
            elif len(func[i]) == 2:
                func_array.append(func[i][0])
    for i in range(len(func)):
        if len(func[i]) == 2:
            func_array.append(func[i][1])
    func_array = sp.Matrix(func_array)
    return func_array

def Jacobian(bus, func_array, sym_variables):
    J = np.zeros([len(func_array), len(func_array)], dtype=object)
    # print(J[0])
    for i in range(len(func_array)):
        for j in range(len(sym_variables)):
            J[i, j] = func_array[i].diff(sym_variables[j])
    J = sp.Matrix(J)
    return J

def PQ_given(bus):
    P_given = np.zeros(len(bus))
    Q_given = np.zeros(len(bus))
    for i in range(len(bus)):
        if bus['Type'][i] == 'PQ':
            P_given[i] = -bus['PL'][i]
            Q_given[i] = -bus['QL'][i]
        elif bus['Type'][i] == 'PV':
            P_given[i] = bus['PG'][i]

    U_given = []
    for i in range(len(bus)):
        if bus['Type'][i] == 'PQ':
            U_given.append(P_given[i])
        elif bus['Type'][i] == 'PV':
            U_given.append(P_given[i])
    for i in range(len(bus)):
        if bus['Type'][i] == 'PQ':
            U_given.append(-Q_given[i]) #확인 필요
    # print(f"U_given : {U_given}")
    return U_given