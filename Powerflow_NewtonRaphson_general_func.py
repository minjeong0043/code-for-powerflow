import pandas as pd
import numpy as np
import sympy as sp
import time
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
                    # print(f"k : {k}")
                    # Y[i, j] += branch['G'][k] + 1j * branch['B'][k]
                    Y[i,j] += 1/(branch['R (pu)'][k] + 1j * branch['X (pu)'][k]) + 1j* branch['B (pu)'][k]/2
                Y_mag[i, j] = np.abs(Y[i, j])
                Y_degrees[i, j] = np.degrees(np.angle(Y[i, j]))
            elif i != j:
                index = np.where(((branch['From'] == i + 1) & (branch['To'] == j + 1)) | ((branch['From'] == j + 1) & (branch['To'] == i + 1)))[0][0]
                # print(f"i, j { i+1, j+1}         index : {index}")
                if index.size == 0:
                    # print("FF")
                    Y[i, j] = 0
                    Y_mag[i, j] = np.abs(Y[i, j])
                    Y_degrees[i, j] = np.degrees(np.angle(Y[i, j]))
                else:
                    Y[i, j] = - 1/(branch['R (pu)'][index] + 1j * branch['X (pu)'][index])
                    Y_mag[i, j] = np.abs(Y[i, j])
                    Y_degrees[i, j] = np.degrees(np.angle(Y[i, j]))
    # print(f"Y_mag : {Y_mag}         Y_degrees: {sp.rad(Y_degrees)}")
    Y_rad = sp.rad(Y_degrees)
    return Y_mag, Y_rad

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
                    # print(f"index: {index}                          index_trans : {index_trans}              set(index)&set(index_trans) : {set(index)&set(index_trans)}")
                    for k in index:
                        Y[i, j] += 1 / (branch['R (pu)'][k] + 1j * branch['X (pu)'][k]) + 1j * branch['B (pu)'][k] / 2
                elif len(set(index) & set(index_trans)) != 0:
                    # print(f"index: {index}                          index_trans : {index_trans}              set(index)&set(index_trans) : {set(index)&set(index_trans)}")
                    # print(set(index) - set(index_trans))
                    for k in (set(index) - set(index_trans)):
                        # print(f"k : {k}")
                        Y[i, j] += 1 / (branch['R (pu)'][k] + 1j * branch['X (pu)'][k]) + 1j * branch['B (pu)'][k] / 2
                    for k in (set(index) & set(index_trans)):
                        # print(f"k & : {k}           branch['B (pu)'][k] : {branch['B (pu)'][k]}")
                        # print(trans)
                        # print([item for item in trans if item[3] == k][0][2])
                        t = [item for item in trans if item[3] == k][0][2]
                        # print(f"t : {t}          np.power(t) : {np.power(t, 2)}")
                        Y[i, j] += (1 / (branch['R (pu)'][k] + 1j * branch['X (pu)'][k])) / np.power(t, 2)
                # for k in index:
                #     Y[i,j] += 1/ (branch['R (pu)'][k] + 1j*branch['X (pu)'][k]) + 1j*branch['B (pu)'][k]/2
            elif i != j:
                index = np.where(((branch['From'] == i + 1) & (branch['To'] == j + 1)) | (
                            (branch['From'] == j + 1) & (branch['To'] == i + 1)))[0]
                # print(f"i+1,j+1: {i+1, j+1}            index: {index}                          index_trans : {index_trans}              set(index)&set(index_trans) : {set(index) & set(index_trans)}")
                if len(set(index) & set(index_trans)) == 0:
                    # print(f"index: {index}                          index_trans : {index_trans}              set(index)&set(index_trans) : {set(index)&set(index_trans)}")
                    if index.size == 0:
                        Y[i, j] = 0
                    else:
                        Y[i, j] = -1 / (branch['R (pu)'][index[0]]) + 1j * branch['X (pu)'][index[0]]
                elif len(set(index) & set(index_trans)) != 0:
                    # print(f"index: {index}                          index_trans : {index_trans}              set(index)&set(index_trans) : {set(index) & set(index_trans)}")
                    for k in (set(index) - set(index_trans)):
                        # print(f"k : {k}")
                        Y[i, j] = -1 / (branch['R (pu)'][k]) + 1j * branch['X (pu)'][k]
                    for k in (set(index) & set(index_trans)):
                        # print(f"k : {k}")
                        t = [item for item in trans if item[3] == k][0][2]
                        # print(f"t : {t}          np.power(t) : {np.power(t, 2)}")
                        # print(trans)
                        # print(t)
                        Y[i, j] = (-1 / (branch['R (pu)'][k] + 1j * branch['X (pu)'][k])) / t
    # print(Y)
    for i in range(len(bus)):
        for j in range(len(bus)):
            Y_mag[i, j] = np.abs(Y[i, j])
            Y_degrees[i, j] = np.degrees(np.angle(Y[i, j]))

    Y_rad = sp.rad(Y_degrees)
    return Y_mag, Y_rad

def DefineVariable(bus_pu, generator_pu):
    V_mag_init = bus_pu['Vm (pu)'].copy()
    delta_init = bus_pu['Va (degree)'].copy()
    for i in range(len(generator_pu)):
        V_mag_init.loc[generator_pu['Bus'][i] - 1] = generator_pu['Voltage setpoint (pu)'][i]
    # print(V_mag_init)

    variables = []
    indices_delta = []
    for i in range(len(bus_pu)):
        if bus_pu['Type'][i] == 'PQ' or bus_pu['Type'][i] == 'PV':
            variables.append(f"delta{i}")
            indices_delta.append(i)
    for i in range(len(bus_pu)):
        if bus_pu['Type'][i] == 'PQ':
            variables.append(f"V{i}")
    sym_variables = sp.symbols(variables)
    # print(f"sym_variables = {sym_variables}")

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

def PQ_given(bus_pu, generator_pu):
    PG = np.zeros(len(bus_pu))
    for i in range(len(generator_pu)):
        PG[generator_pu['Bus'][i] -1] = generator_pu['PG (pu)'][i]
        # print(f"PG : {PG}")

    P_given = np.empty(len(bus_pu))
    Q_given = np.empty(len(bus_pu))
    for i in range(len(bus_pu)):
        if bus_pu['Type'][i] == 'PQ':
            P_given[i] = -bus_pu['Pload (pu)'][i]
            Q_given[i] = -bus_pu['Qload (pu)'][i]
        elif bus_pu['Type'][i] == 'PV':
            P_given[i] = PG[i] - bus_pu['Pload (pu)'][i]
            # print(f" i : {i},   PG {PG},       bus_pu :{isinstance(bus_pu['Pload (pu)'][i], np.float64())}")

    U_given = []
    for i in range(len(bus_pu)):
        if bus_pu['Type'][i] == 'PQ':
            U_given.append(P_given[i])
        elif bus_pu['Type'][i] == 'PV':
            U_given.append(P_given[i])
    for i in range(len(bus_pu)):
        if bus_pu['Type'][i] == 'PQ':
            U_given.append(Q_given[i])
    return U_given

def NewtonRaphson(x, Ufunc, bus, U_given, J, indices_delta):
    Ufunc_cal = Ufunc.subs(x)
    # print(f"Ufunc_cal :{Ufunc_cal},    U_given : {U_given}")
    del_U = np.zeros(len(U_given))
    for i in range(len(U_given)):
        del_U[i] = U_given[i] - Ufunc_cal[i]
    del_U = sp.Matrix(del_U)
    # print(f"de_U : {del_U}")
    J_val = J.subs(x)
    # print(f"matrix : {J_val @ del_U}")
    # print(f"J_val.inv() : {J_val.inv()}")
    del_x = J_val.inv() @ del_U
    # print(f"del_x : {del_x}")
    del_x_degree = del_x.copy()
    for i in range(len(indices_delta)):
        del_x_degree[i] = del_x[i] * 180 /np.pi

    a = enumerate(x)
    x_new = {}
    x_new_degree = {}
    x = x.values()
    x_values = list(x)
    # for i, key in enumerate(var0_dict):
    for i, key in a:
        x_new[key] = np.array(x_values[i]) + del_x[i]
        if 'delta' in str(key):
            x_new_degree[key] = (np.array(x_values[i]) + del_x[i]) * 180 / np.pi
        else:
            x_new_degree[key] = np.array(x_values[i]) + del_x[i]

    return x_new, x_new_degree

def NewtonRaphson2(x, Ufunc, bus, U_given, J, indices_delta): # do not have J.inv()
    start_time = time.time()
    Ufunc_cal = Ufunc.subs(x)
    del_U = np.zeros(len(U_given))
    for i in range(len(U_given)):
        del_U[i] = U_given[i] - Ufunc_cal[i]
    del_U = sp.Matrix(del_U)
    # print(f"del_U :{del_U}")
    J_val = J.subs(x)
    print(f"J_val : {J_val}")
    # GaussElimination
    # J_val = J_val.applyfunc(lambda o: sp.Rational(o))
    del_U_updated = sp.zeros(len(del_U), 1)
    J_val_updated = sp.zeros(len(del_U), len(del_U))
    J_val_updated[0, :] = J_val[0, :]
    del_U_updated[0] = del_U[0]
    for i in range(len(del_U)):
        for j in range(i + 1, len(del_U)):
            k = J_val[i, j] / J_val[i, i]
            J_val_updated[j, :] = J_val.row(j) - k * J_val.row(i)
            del_U_updated[j] = del_U[j] - k * del_U[i]
    # J_val_updated = J_val_updated.evalf()
    print(f"J_val+_updated :{J_val_updated}")
    print(f"del_U_updated  :{del_U_updated}")

    del_x = sp.zeros(len(del_U), 1)
    del_x[-1] = 1 / J_val_updated[-1, -1] * del_U[-1]
    print(f"-------------")
    print(" -- %s seconds --- " % (time.time() - start_time))

    for i in range(len(del_U) - 2, -1, -1):
        p = 0
        # for j in range(i + 1, len(del_U)):
            # p += J_val_updated[i, j] * del_x[j]
        p = np.dot(J_val_updated[i,i+1:], del_x[i+1:])
        del_x[i] = 1 / J_val_updated[i, i] * (del_U_updated[i] - p)
    print(f"del_x : {del_x}")
    # del_x_1 = J_val.inv() @ del_U
    # print(f"del_x :{del_x},          del_x_1 :{del_x_1}")
    print(f"-------------")
    print(" -- %s seconds --- " % (time.time() - start_time))
    del_x_degree = del_x.copy()
    for i in range(len(indices_delta)):
        del_x_degree[i] = del_x[i] * 180 /np.pi

    a = enumerate(x)
    x_new = {}
    x_new_degree = {}
    x = x.values()
    x_values = list(x)
    # for i, key in enumerate(var0_dict):
    for i, key in a:
        x_new[key] = np.array(x_values[i]) + del_x[i]
        if 'delta' in str(key):
            x_new_degree[key] = (np.array(x_values[i]) + del_x[i]) * 180 / np.pi
        else:
            x_new_degree[key] = np.array(x_values[i]) + del_x[i]
    print(f"x_new : {x_new}")
    return x_new, x_new_degree


def NewtonRaphson3(x, Ufunc, bus, U_given, J, indices_delta):
    # start_time = time.time()
    Ufunc_cal = Ufunc.subs(x).evalf()
    del_U = np.zeros(len(U_given))

    # Ufunc와 U_given 차이를 계산
    for i in range(len(U_given)):
        del_U[i] = float(U_given[i] - Ufunc_cal[i])
    del_U = sp.Matrix(del_U)

    # Jacobi 행렬 계산
    J_val = J.subs(x).evalf()
    J_val = np.array(J_val).astype(np.float64)  # Sympy 행렬을 numpy 배열로 변환
    del_U = np.array(del_U).astype(np.float64)  # Sympy 행렬을 numpy 배열로 변환
    print(f"                    J_val : {J_val} {type(J_val)}\ndel_U :{del_U} {type(del_U)}")
    # LU 분해 적용 (scipy 라이브러리 사용)
    lu, piv = lu_factor(J_val)  # LU 분해
    del_x = lu_solve((lu, piv), del_U)  # LU 분해 결과로 연립 방정식 풀이

    # print(f"del_x : {del_x}")
    # print(" -- %s seconds --- " % (time.time() - start_time))

    # 각 변수의 변화량을 라디안에서 각도로 변환
    # del_x_degree = del_x.copy()
    # for i in range(len(indices_delta)):
    #     del_x_degree[i] = del_x[i] * 180 / np.pi

    # 새로운 값 계산
    a = enumerate(x)
    x_new = {}
    # x_new_degree = {}
    x_values = list(x.values())
    print(f"del_x[0] : {del_x[0]} {type(del_x[0])}")
    # 각 변수에 대해 업데이트된 값을 계산
    for i, key in a:
        x_new[key] = np.array(x_values[i]) + del_x[i][0]
        # print(f"np.array(x_values[i]) : {np.array(x_values[i])} {type(np.array(x_values[i]))}, {x_values[i]}, del_x : {del_x[i][0]}{type(del_x[i])}")
        # if 'delta' in str(key):
        #     x_new_degree[key] = (np.array(x_values[i]) + del_x[i][0]) * 180 / np.pi
        # else:
        #     x_new_degree[key] = np.array(x_values[i]) + del_x[i][0]

    print(f"x_new : {x_new}")
    return x_new, del_x

def Result_values(x, bus_pu, V_mag, delta, Y_mag, Y_rad):
    V_result = np.empty(len(bus_pu), dtype=float)
    delta_result = np.empty(len(bus_pu), dtype=float)
    for i in range(len(bus_pu)):
        if isinstance(V_mag[i], np.float64):
            V_result[i] = V_mag[i]
        else:
            V_result[i] = x.get(V_mag[i], np.nan)
        if isinstance(delta[i], np.float64) or isinstance(delta[i], np.int64):
            delta_result[i] = delta[i]
        else:
            delta_result[i] = x.get(delta[i], np.nan)

    P_result = np.zeros(len(bus_pu), dtype=float)
    Q_result = np.zeros(len(bus_pu), dtype=float)

    # P, Q 계산
    for i in range(len(bus_pu)):
        if bus_pu['Type'][i] == 'PQ':
            for k in range(len(bus_pu)):
                angle_diff = float(delta_result[i] - delta_result[k] - Y_rad[i, k])
                P_result[i] += Y_mag[i, k] * V_result[k] * np.cos(angle_diff)
                Q_result[i] += Y_mag[i, k] * V_result[k] * np.sin(angle_diff)
            P_result[i] = V_result[i] * P_result[i]
            Q_result[i] = V_result[i] * Q_result[i]

        elif bus_pu['Type'][i] == 'PV':
            for k in range(len(bus_pu)):
                angle_diff = float(delta_result[i] - delta_result[k] - Y_rad[i, k])
                P_result[i] += Y_mag[i, k] * V_result[k] * np.cos(angle_diff)
                Q_result[i] += Y_mag[i, k] * V_result[k] * np.sin(angle_diff)
            P_result[i] = V_result[i] * P_result[i]
            Q_result[i] = V_result[i] * Q_result[i]

        elif bus_pu['Type'][i] == 'Swing':
            for k in range(len(bus_pu)):
                angle_diff = float(delta_result[i] - delta_result[k] - Y_rad[i, k])
                P_result[i] += Y_mag[i, k] * V_result[k] * np.cos(angle_diff)
                Q_result[i] += Y_mag[i, k] * V_result[k] * np.sin(angle_diff)
            P_result[i] = V_result[i] * P_result[i]
            Q_result[i] = V_result[i] * Q_result[i]

    # 결과 출력
    print(f"V_result     : {V_result}")
    print(f"delta_result : {delta_result}")
    print(f"P_result     : {P_result}")
    print(f"Q_result     : {Q_result}")
    return V_result, delta_result, P_result, Q_result