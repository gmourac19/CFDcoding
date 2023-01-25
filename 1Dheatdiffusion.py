import numpy as np
import matplotlib.pyplot as mp

A = 160*10**(-4)
k = 20
L = 0.6*10**(-2)
deltax = 0.05*10**(-2)
SV = 0
n = int(L / deltax)  # no. of centroids
n_aux = n - 1
bc_left = -800  # W
T_bc2 = 85  # C

# matrix CT = b
C = np.zeros([n, n])
b = np.zeros([n, 1])
# print(C)
# print(b)

# ap, al, ar for m = 1
al_left = 0
ar_left = k * A / deltax
Su_left = -bc_left
ap_left = al_left + ar_left
C[0, 0] = ap_left
C[0, 1] = -ar_left
b[0] = Su_left
# print(C)

# ap, al, ar for 1<m<n
al_mid = k*A/deltax
ar_mid = k*A/deltax
ap_mid = al_mid + ar_mid
for index in range(1, n_aux):
    C[index, index-1] = -al_mid
    C[index, index] = ap_mid
    C[index, index+1] = -ar_mid
    b[index] = 100  # geração interna de energia

# ap, al, ar for m = n
al_n = k*A/deltax
ar_n = 0
ap_n = al_n + 2*k*A/deltax
Su_n = T_bc2*(2*k*A/deltax)
C[n_aux, n_aux-1] = -al_n
C[n_aux, n_aux] = ap_n
b[n_aux] = Su_n


res = 0
i = 0
j = 0
T = np.zeros([n, 1])
c = np.zeros([n, 1])
T1 = np.zeros([n, n])
T2 = np.zeros([n, n])

# Jacobi method as presented in Malalasekera & Versteeg
# for iteration in range(0, 10000):
#     T_ant = T
#     for i in range(0, n):
#         for j in range(0, n):
#             if i != j:
#                 T1[i, j] = -C[i, j] / C[i, i]
#             else:
#                 T1[i, j] = 0
#             c[i] = b[i] / C[i, i]
#
#     T = np.dot(np.array(T1), np.array(T_ant)) + np.array(c)
#     res = abs(np.dot(np.array(np.transpose(T)), np.array(T)) - np.dot(np.array(np.transpose(T_ant)), np.array(T_ant)))
#     print('residuals (iteration {}): {}'.format(iteration, res))
#     if res > 1e-5:
#         continue
#     else:
#         break

#Gauss-Seidel method as presented in Malalasekera & Versteeg
for iteration in range(0, 10000):
    T_ant = T
    for i in range(0, n):
        for j in range(0, n):
            if i > j:
                T1[i, j] = -C[i, j] / C[i, i]
                T2[i, j] = 0
            elif i == j:
                T1[i, j] = 0
                T2[i, j] = 0
            elif i < j:
                T1[i, j] = 0
                T2[i, j] = -C[i, j] / C[i, i]
            c[i] = b[i] / C[i, i]

    T = np.dot(np.array(T1), np.array(T)) + np.dot(np.array(T2), np.array(T_ant)) + np.array(c)
    res = abs(np.dot(np.array(np.transpose(T)), np.array(T)) - np.dot(np.array(np.transpose(T_ant)), np.array(T_ant)))
    print('residuals (iteration #{}): {}'.format(iteration, res))
    if res > 1e-5:
        continue
    else:
        break

# print(T1)
# print(c)
# print(T)
# print(res)
# T = np.linalg.solve(C,b)
# print(C)
x = 0
posx = []
for ind in range(0, n):
    posx.append((ind+1)*deltax-deltax)
print((posx))
mp.plot(posx, T, marker='.')
mp.show()
