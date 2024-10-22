import pandas as pd
import numpy as np
import cmath

data_name = 'ac_case25.xlsx'

# read data
bus = pd.read_excel(data_name, sheet_name='bus')
branch = pd.read_excel(data_name, sheet_name='branch')
transformer = pd.read_excel(data_name, sheet_name='transformer')
generator = pd.read_excel(data_name, sheet_name='generator')
Sbase = pd.read_excel(data_name, sheet_name='param', header=None)[1][0]

# check Swing bus
index_Swing = np.where(bus['Type'] == 'Swing')[0]

# Change Unit from MVAR, MW to pu
bus_pu = bus.copy()
bus_pu['Pload (MW)'] = bus['Pload (MW)'] / Sbase
bus_pu['Qload (MVAR)'] = bus['Qload (MVAR)'] / Sbase
bus_pu = bus_pu.rename(columns = {'Pload (MW)': 'Pload (pu)', 'Qload (MVAR)' : 'Qload (pu)'})

# Cal Ybus
# Ybus considering Transformer
trans = []
index_trans = []
for i in range(len(transformer)):
    for j in range(len(branch)):
        if (((branch['From'][j] == transformer['From'][i]) & (branch['To'][j] == transformer['To'][i])) | ((branch['From'][j] == transformer['To'][i]) & (branch['To'][j] == transformer['From'][i]))):
            # From, To, tap, index in branch sheet
            trans.append([transformer['From'][i], transformer['To'][i], transformer['Tap'][i], j])
            index_trans.append(j)

Y = np.zeros([len(bus), len(bus)], dtype=complex)
for i in range(len(bus)):
    for j in range(len(bus)):
        if i == j:
            index = np.where((branch['From'] == i+1) | (branch['To'] == i+1))[0]
            if len(set(index)&set(index_trans)) == 0:
                # print(f"index: {index}                          index_trans : {index_trans}              set(index)&set(index_trans) : {set(index)&set(index_trans)}")
                for k in index:
                    Y[i, j] += 1 / (branch['R (pu)'][k] + 1j * branch['X (pu)'][k]) + 1j * branch['B (pu)'][k] / 2
            elif len(set(index)&set(index_trans)) != 0:
                # print(f"index: {index}                          index_trans : {index_trans}              set(index)&set(index_trans) : {set(index)&set(index_trans)}")
                # print(set(index) - set(index_trans))
                for k in (set(index) - set(index_trans)):
                    # print(f"k : {k}")
                    Y[i, j] += 1 / (branch['R (pu)'][k] + 1j * branch['X (pu)'][k]) + 1j * branch['B (pu)'][k] / 2
                for k in  (set(index) & set(index_trans)):
                    # print(f"k & : {k}           branch['B (pu)'][k] : {branch['B (pu)'][k]}")
                    # print(trans)
                    # print([item for item in trans if item[3] == k][0][2])
                    t = [item for item in trans if item[3] == k][0][2]
                    # print(f"t : {t}          np.power(t) : {np.power(t, 2)}")
                    Y[i,j] += (1/ (branch['R (pu)'][k] + 1j*branch['X (pu)'][k])) / np.power(t, 2)
            # for k in index:
            #     Y[i,j] += 1/ (branch['R (pu)'][k] + 1j*branch['X (pu)'][k]) + 1j*branch['B (pu)'][k]/2
        elif i != j:
            index = np.where(((branch['From'] == i+1) & (branch['To'] == j+1)) | ((branch['From'] == j+1) & (branch['To'] == i+1)))[0]
            # print(f"i+1,j+1: {i+1, j+1}            index: {index}                          index_trans : {index_trans}              set(index)&set(index_trans) : {set(index) & set(index_trans)}")
            if len(set(index) & set(index_trans)) == 0:
                # print(f"index: {index}                          index_trans : {index_trans}              set(index)&set(index_trans) : {set(index)&set(index_trans)}")
                if index.size == 0:
                    Y[i,j] = 0
                else:
                    Y[i,j] = -1/ (branch['R (pu)'][index[0]]) + 1j* branch['X (pu)'][index[0]]
            elif len(set(index) & set(index_trans)) != 0:
                # print(f"index: {index}                          index_trans : {index_trans}              set(index)&set(index_trans) : {set(index) & set(index_trans)}")
                for k in (set(index) - set(index_trans)):
                    print(f"k : {k}")
                    Y[i,j] = -1/ (branch['R (pu)'][k]) + 1j* branch['X (pu)'][k]
                for k in  (set(index) & set(index_trans)):
                    # print(f"k : {k}")
                    t = [item for item in trans if item[3] == k][0][2]
                    # print(f"t : {t}          np.power(t) : {np.power(t, 2)}")
                    # print(trans)
                    # print(t)
                    Y[i,j] = (-1/ (branch['R (pu)'][k] + 1j* branch['X (pu)'][k])) / t

print(Y)
