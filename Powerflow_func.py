import numpy as np
import pandas as pd

def ReadData_Unit_pu(data_path):
    # read input file and define parameters
    bus = pd.read_excel(data_path, sheet_name='bus')
    branch = pd.read_excel(data_path, sheet_name='branch')
    transformer = pd.read_excel(data_path, sheet_name='transformer')
    generator = pd.read_excel(data_path, sheet_name='generator')
    Sbase = pd.read_excel(data_path, sheet_name='param', header=None).iloc[0, 1]

    # Transform Unit to pu in bus sheet
    bus['Pload (MW)'] = bus['Pload (MW)'] / Sbase
    bus['Qload (MVAR)'] = bus['Qload (MVAR)'] / Sbase
    bus_pu = bus.rename(columns={'Pload (MW)': 'Pload (pu)', 'Qload (MVAR)': 'Qload (pu)'})
    # Transform Unit to pu in generator sheet
    generator['PG (MW)'] = generator['PG (MW)'] / generator['MBASE (MW)']
    generator['QG (MVAR)'] = generator['QG (MVAR)'] / generator['MBASE (MW)']
    generator['QMAX (MVAR)'] = generator['QMAX (MVAR)'] / generator['MBASE (MW)']
    generator['QMIN (MVAR)'] = generator['QMIN (MVAR)'] / generator['MBASE (MW)']
    generator_pu = generator.rename(columns={'PG (MW)': 'PG (pu)', 'QG (MVAR)': 'QG (pu)', 'QMAX (MVAR)': 'QMAX (pu)', 'QMIN (MVAR)': 'QMIN (pu)'})

    return bus_pu, branch, transformer, generator_pu, Sbase

def TransformUnit(bus, Sbase):
    bus['Pload (MW)'] = bus['Pload (MW)'] / Sbase
    bus['Qload (MVAR)'] = bus['Qload (MVAR)'] / Sbase
    bus_pu = bus.rename(columns = {'Pload (MW)' : 'Pload (pu)', 'Qload (MVAR)' : 'Qlaod (pu)'})
    return bus_pu

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

def Ybus_considering_trans(bus, branch, transformer):
    size_Mtx = len(bus)
    Y = np.zeros([size_Mtx, size_Mtx], dtype=complex)

    trans = []
    index_trans = []
    tap = np.zeros(len(branch))
    for i in range(len(transformer)):
        trans.append([transformer['From'][i], transformer['To'][i], transformer['Tap'][i]])

    # 변압기 부분 branch에서의 index와 tap 재배치
    for i in range(len(transformer)):
        a = np.where((branch['From'] == trans[i][0]) & (branch['To'] == trans[i][1]))[0]
        b = np.where((branch['From'] == trans[i][1]) & (branch['To'] == trans[i][0]))[0]
        s = set(a) | set(b)
        for k in s:
            tap[k] = trans[i][2]
        if len(s) != 0:
            index_trans.append(s)
    index_trans = [next(iter(s)) for s in index_trans] # set to list

    #어드미턴스 계싼
    for i in range(size_Mtx):
        for j in range(size_Mtx):
            if i == j:
                index = np.where((branch['From'] == i+1) | (branch['To'] == i+1))[0]
                index_without_trans = set(index) - set(index_trans)
                index_with_trans = set(index) & set(index_trans)

                if len(index_without_trans) == 0:
                    Y[i,j] += 0
                else:
                    for k in index_without_trans:
                        Y[i,j] += 1/(branch['R (pu)'][k] + 1j*branch['X (pu)'][k]) + 1j*branch['B (pu)'][k]/2

                if len(index_with_trans) == 0:
                    Y[i,j] += 0
                else:
                    for k in index_with_trans:
                        Y[i,j] += (1/(branch['R (pu)'][k] + 1j*branch['X (pu)'][k])) / np.power(tap[k], 2)
            elif i != j:
                index = np.where(((branch['From'] == i + 1) & (branch['To'] == j + 1)) | ((branch['From'] == j + 1) & (branch['To'] == i + 1)))[0]
                index_without_trans = set(index) - set(index_trans)
                index_with_trans = set(index) & set(index_trans)
                # if len(index_with_trans) > 1 or len(index_without_trans) > 1: # 두 번씩 나타나는 애들
                #     print(f"i, j = {i, j},            index_without_trans = {index_without_trans},                   index_with_trans = {index_with_trans}")
                if (len(index_without_trans) != 0) and (len(index_with_trans) == 0):
                    Y[i, j] = - 1 / (branch['R (pu)'][next(iter(index_without_trans))] + 1j * branch['X (pu)'][
                        next(iter(index_without_trans))])
                elif (len(index_with_trans) != 0) and (len(index_without_trans) == 0):
                    Y[i, j] = - (1 / (branch['R (pu)'][next(iter(index_with_trans))] + 1j * branch['X (pu)'][
                        next(iter(index_with_trans))])) / tap[next(iter(index_with_trans))]
                    # print(f"i, j = {i, j},            index_without_trans = {index_without_trans},                   index_with_trans = {index_with_trans}")
                elif (len(index_with_trans) == 0) and (len(index_without_trans) == 0):
                    # print(f"i, j = {i, j},            index_without_trans = {index_without_trans},                   index_with_trans = {index_with_trans}")
                    Y[i, j] = 0
                else:
                    print("WWWWWWWWWWWWWOOOOOOOOOOOOOOOOOWWWWWWWW")

    return Y

def Ybus_considering_trans2(bus, branch, transformer):
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
                        print(f"k : {k}")
                        Y[i, j] = -1 / (branch['R (pu)'][k]) + 1j * branch['X (pu)'][k]
                    for k in (set(index) & set(index_trans)):
                        # print(f"k : {k}")
                        t = [item for item in trans if item[3] == k][0][2]
                        # print(f"t : {t}          np.power(t) : {np.power(t, 2)}")
                        # print(trans)
                        # print(t)
                        Y[i, j] = (-1 / (branch['R (pu)'][k] + 1j * branch['X (pu)'][k])) / t
    return Y

def Cal_PQ(V, Y, i, size_bus): # 전압 크기, 위상, 어드미턴스, bus i에서의 계산 PQ 전력
    # V_value, delta = cmath.polar(V)
    V_value = []
    delta = []

    # V 배열의 각 요소에 대해 크기와 각도 계산
    # for v in V:
    #     v_val, v_delta = cmath.polar(v)
    #     V_value.append(v_val)
    #     delta.append(v_delta)

    V_value = abs(V)
    delta = np.angle(V)

    P_cal = 0
    Q_cal = 0
    # t = 0
    for j in range(size_bus):
        P_cal += V_value[i] * V_value[j] * (Y[i][j].real * np.cos(delta[i] - delta[j]) + Y[i][j].imag * np.sin(delta[i] - delta[j]))
        Q_cal += V_value[i] * V_value[j] * (Y[i][j].real * np.sin(delta[i] - delta[j]) - Y[i][j].imag * np.cos(delta[i] - delta[j]))

    # print("P_cal : ", P_cal)
    # print("Q_cal : ", Q_cal)
    return P_cal, Q_cal

def GaussSeidel_PQ(data_path, bus_index, V, Y):
    bus_pu, branch, transformer, generator_pu, Sbase = ReadData_Unit_pu(data_path)
    PG = np.zeros(len(bus_pu))
    QG = np.zeros(len(bus_pu))
    for i in range(len(generator_pu)):
        # print(f"generator_pu['Bus'][i] : {generator_pu['Bus'][i]}, PG : {generator_pu['PG (pu)'][i]}")
        # PG[generator_pu['Bus'][i] - 1] += generator_pu['PG (pu)'][i]
        # QG[generator_pu['Bus'][i] - 1] += generator_pu['QG (pu)'][i]
        PG[generator_pu['Bus'][i] - 1] = generator_pu['PG (pu)'][i]
        QG[generator_pu['Bus'][i] - 1] = generator_pu['QG (pu)'][i]

    # print(f"PG : {PG}")


    P_given = PG - bus_pu['Pload (pu)']
    Q_given = QG - bus_pu['Qload (pu)']

    tolerance = 1e-6
    max_iter = 100
    iteration = 0
    converged = False

    while not converged and iteration < max_iter:
        P_cal, Q_cal = Cal_PQ(V, Y, bus_index, len(bus_pu))
        del_P = abs(P_given[bus_index] - P_cal)
        del_Q = abs(Q_given[bus_index] - Q_cal)
        # print(f"del_P : {del_P},                          del_Q : {del_Q},                 V : {V[bus_index]}")
        if ((del_P < tolerance) and (del_Q < tolerance)):
            converged = True
            # print(f"Converged at i = {iteration}, V = {V[bus_index]}")
            break

        b = 0
        for j in range(len(bus_pu)):
            if bus_index != j: # 자기 자신이 아닌 거에서만
                b += Y[bus_index,j] * V[j]
        # V[bus_index] = (1/Y[bus_index, bus_index]) * (((P_given[bus_index] - 1j* Q_given[bus_index]) / V[bus_index].conjugate()) - b)
        V_new = (1/Y[bus_index, bus_index]) * (((P_given[bus_index] - 1j* Q_given[bus_index]) / V[bus_index].conjugate()) - b)
        V_mag = np.abs(V_new)
        if V_mag < bus_pu['minVm'][bus_index]:
            print(f"Voltage magnitude is out of range : {V_mag}")
            V[bus_index] = bus_pu['minVm'][bus_index] * np.exp(1j * np.angle(V_new)) # 이렇게 해도 되나 생각해보자
        elif V_mag > bus_pu['maxVm'][bus_index]:
            print(f"Voltage magnitude is out of range : {V_mag}")
            V[bus_index] = bus_pu['maxVm'][bus_index] * np.exp(1j * np.angle(V_new))
        else:
            V[bus_index] = V_new

        iteration += 1
        # print(f"iteration = {iteration}")

    return V

def GaussSeidel_PV(data_path, bus_index, V, Y):
    bus_pu, branch, transformer, generator_pu, _ = ReadData_Unit_pu(data_path)
    PG = np.zeros(len(bus_pu))
    QG = np.zeros(len(bus_pu))
    for i in range(len(generator_pu)):
        # print(f"generator_pu['Bus'][i] : {generator_pu['Bus'][i]}")
        # PG[generator_pu['Bus'][i] - 1] += generator_pu['PG (pu)'][i]
        # QG[generator_pu['Bus'][i] - 1] += generator_pu['QG (pu)'][i]
        PG[generator_pu['Bus'][i] - 1] = generator_pu['PG (pu)'][i]
        QG[generator_pu['Bus'][i] - 1] = generator_pu['QG (pu)'][i]

    P_given = PG - bus_pu['Pload (pu)']
    Q = 0
    V_given = abs(V)

    tolerance = 1e-6
    max_iter = 100
    iteration = 0
    converged = False

    while not converged and iteration < max_iter:
        P_cal, Q_cal = Cal_PQ(V, Y, bus_index, len(bus_pu))
        del_P = abs(P_given[bus_index] - P_cal)
        # del_V = abs(abs(V_given) - abs(V[bus_index]))
        # del_V = abs(V_given[bus_index] - abs(V[bus_index]))
        del_V = abs(V_given[bus_index] - abs(V[bus_index]))
        Q = Q_cal
        # print(f"del_P : {del_P},                          del_V : {del_V},                 V : {V[bus_index]}")
        if ((del_P < tolerance) and (del_V < tolerance)):
            converged = True
            # print(f"Converged at i = {iteration}, V = {V[bus_index]}")
            break

        b = 0
        for j in range(len(bus_pu)):
            if bus_index != j:
                b += Y[bus_index,j] * V[j]
        V_new = (1/Y[bus_index, bus_index]) * (((P_given[bus_index] - 1j* Q) / V[bus_index].conjugate()) - b)
        #위 식 다른 거로도 돌려보고 결과 비교해보기
        # print("!")
        # print(f"type : {type(bus_index)}        {type(V_given)}            {type(V_new)}         {type(V)}")
        # print(f'           {V_given}              {V}')
        # V_given = np.array(V_given)
        # V = np.array(V)
        V[bus_index] = V_given[bus_index] * np.exp(1j*np.angle(V_new))
        # print("1!!!!!!!!!!!!")
        iteration += 1
        # print(f"iteration = {iteration}")

    return V


