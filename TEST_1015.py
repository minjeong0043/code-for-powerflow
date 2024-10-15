import numpy as np
import pandas as pd

# 입력 데이터 로드
bus = pd.read_excel("ac_case25_example.xlsx", sheet_name='bus')
branch = pd.read_excel("ac_case25_example.xlsx", sheet_name='branch')

# print(branch)
# 행렬 크기 결정 (bus의 개수)
# size_Mtx = len(bus)
size_Mtx = 5
# Ybus 행렬 초기화 (복소수 행렬)
Y = np.zeros([size_Mtx, size_Mtx], dtype=complex)

# branch 데이터에서 필요한 값 추출
From = branch['From'] - 1  # 버스 번호가 1부터 시작할 경우 0부터 시작하도록 맞춤
To = branch['To'] - 1  # 동일하게 0부터 시작하도록 맞춤
R = branch['R (pu)']
X = branch['X (pu)']
B = branch['B (pu)']

for i in range(size_Mtx):
    for j in range(size_Mtx):
        if i == j:
            index = np.where(From == i)[0]
            # print("index : ", index)
            if index.size == 0:
                Y[i, j] = 0
            else:
                for k in index:
                    Y[i, j] += 1/(R[k] + 1j * X[k]) + 1j*B[k]/2
        else:
            index = np.where(((From == i) & (To == j))| ((From == j) & (To == i)))[0]
            if index.size == 0:
                Y[i, j] = 0
            else:
                Y[i, j] = - 1/ (R[index] + 1j *X[index])
                # for k in index:
                #     Y[i,j] -= 1/(R[k] + 1j*X[k])

print(Y)

# Ybus 계산
# for i in range(len(branch)):
#     # 어드미턴스 계산 (Y = 1 / (R + jX))
#     y = 1 / (R[i] + 1j * X[i])
#
#     # Shunt 용량이 있는 경우
#     b_shunt = 1j * B[i] / 2
#
#     # 대각 성분 (Y[i, i])
#     Y[From[i], From[i]] += y + b_shunt
#     Y[To[i], To[i]] += y + b_shunt
#
#     # 비대각 성분 (Y[i, j]와 Y[j, i], 대칭)
#     Y[From[i], To[i]] -= y
#     Y[To[i], From[i]] -= y
#
# # 결과 출력
# print("Ybus 행렬:")
# print(Y)
# print(Y[1])

# print(Y[1][1])
# y = 1/ (R[0] + 1j*X[0])
# print(y)
V = np.ones(size_Mtx, dtype = complex)
P = np.zeros(size_Mtx, dtype = int)
Q = np.zeros(size_Mtx)

# Data
P[1] = -8 #pu
Q[1] = -2.8 #pu
# 부하 모선 가우스자이델법
# iteration 1
V[1] = (1/Y[1][1]) * (((P[1]-1j*Q[1])/V[0].conjugate()) - (Y[1][3]*V[3] + Y[1][4]*V[4]))
# print(Y[1][1])
print(V[1])
# print("-------------------------",((P[1]-1j*Q[1])/V[0].conjugate()) - (Y[1][3]*V[3] + Y[1][4]*V[4]))
# iteration 2
V[1] = (1/Y[1][1]) * (((P[1]-1j*Q[1])/V[1].conjugate()) - (Y[1][3]*V[3] + Y[1][4]*V[4]))
print(V[1])
