import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

def solve_explicit(rho, Lambda, beta, lambda_, h):
    t = np.arange(0, 1+h, h)
    n = np.zeros(len(t))
    C = np.zeros(len(t))

    # 定义方程组
    def Fun1(t, n, C):
        return (rho - beta) * n / Lambda + lambda_ * C

    def Fun2(t, n, C):
        return beta * n / Lambda - lambda_ * C

    # 定义初值
    n[0] = 1
    C[0] = 1

    # 更新下一个dn/dt和dC/dt
    for i in range(1, len(t)):
        n[i] = n[i-1] + h * Fun1(t[i-1], n[i-1], C[i-1])
        C[i] = C[i-1] + h * Fun2(t[i-1], n[i-1], C[i-1])

    return t, n, C

# 定义参数
rho = 0.0022
Lambda = 10**(-3)
beta = 0.0065
lambda_ = 0.078
h = 0.001

# 调用求解函数
t, n, C = solve_explicit(rho, Lambda, beta, lambda_, h)

# 绘图
plt.figure(figsize=(10, 6))
plt.grid(True)
plt.plot(t, n, 'rh', markersize=3, label='$n(t)$')
plt.plot(t, C, 'bo', markersize=3, label='$C(t)$')
plt.xlabel('t / s')
plt.legend()
plt.title(f'显式法求解中子动力学方程，步长：{h}')
plt.show()