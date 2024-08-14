import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

def solve_RK4(rho, Lambda, beta, lambda_, h):
    t = np.arange(0, 1, h)
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
        K11 = Fun1(t[i-1], n[i-1], C[i-1])
        K21 = Fun2(t[i-1], n[i-1], C[i-1])
        K12 = Fun1(t[i-1] + 0.5*h, n[i-1] + 0.5*h*K11, C[i-1] + 0.5*h*K21)
        K22 = Fun2(t[i-1] + 0.5*h, n[i-1] + 0.5*h*K11, C[i-1] + 0.5*h*K21)
        K13 = Fun1(t[i-1] + 0.5*h, n[i-1] + 0.5*h*K12, C[i-1] + 0.5*h*K22)
        K23 = Fun2(t[i-1] + 0.5*h, n[i-1] + 0.5*h*K12, C[i-1] + 0.5*h*K22)
        K14 = Fun1(t[i-1] + h, n[i-1] + h*K13, C[i-1] + h*K23)
        K24 = Fun2(t[i-1] + h, n[i-1] + h*K13, C[i-1] + h*K23)
        n[i] = n[i-1] + h/6 * (K11 + 2*K12 + 2*K13 + K14)
        C[i] = C[i-1] + h/6 * (K21 + 2*K22 + 2*K23 + K24)

    return t, n, C

# 定义参数
rho = 0.0022
Lambda = 10**(-3)
beta = 0.0065
lambda_ = 0.078
h = 0.0001

# 调用求解函数
t, n, C = solve_RK4(rho, Lambda, beta, lambda_, h)

# 绘图
plt.figure(figsize=(10, 6))
plt.grid(True)
font = FontProperties(fname=r"c:\windows\fonts\simsun.ttc", size=14)  # 引入中文字体
plt.plot(t, n, 'rh', markersize=3, label='$n(t)$')
plt.plot(t, C, 'bo', markersize=3, label='$C(t)$')
plt.xlabel('t / s', fontproperties=font)
plt.legend(prop=font)
plt.title(f'四阶R-K法求解中子动力学方程，步长：{h}', fontproperties=font)
plt.show()