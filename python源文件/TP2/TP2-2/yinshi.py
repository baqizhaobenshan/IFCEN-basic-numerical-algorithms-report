import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

def solve_implicit(rho, Lambda, beta, lambda_, h):
    t = np.arange(0, 1+h, h)
    n = np.zeros(len(t))
    C = np.zeros(len(t))

    # 定义初值
    n[0] = 1
    C[0] = 1

    # 更新下一个dn/dt和dC/dt
    for i in range(1, len(t)):
        A = np.array([[1 - h * (rho - beta) / Lambda, -h * lambda_], [-h * beta / Lambda, 1 + h * lambda_]])
        B = np.array([n[i-1], C[i-1]])
        X = np.linalg.solve(A, B)
        n[i] = X[0]
        C[i] = X[1]

    return t, n, C

# 定义参数
rho = 0.0022
Lambda = 10**(-3)
beta = 0.0065
lambda_ = 0.078
h = 0.01

# 调用求解函数
t, n, C = solve_implicit(rho, Lambda, beta, lambda_, h)

# 绘图
plt.figure(figsize=(10, 6))
plt.grid(True)
font = FontProperties(fname=r"c:\windows\fonts\simsun.ttc", size=14)  # 引入中文字体
plt.plot(t, n, 'rh', markersize=3, label='$n(t)$')
plt.plot(t, C, 'bo', markersize=3, label='$C(t)$')
plt.xlabel('t / s', fontproperties=font)
plt.legend(prop=font)
plt.title('隐式法求解中子动力学方程', fontproperties=font)
plt.show()