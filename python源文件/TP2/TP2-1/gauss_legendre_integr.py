from scipy.integrate import quad
import numpy as np

def gauss_legendre_integr(a, b):
    # 高斯点
    GaussP = [-np.sqrt(15)/5, 0, np.sqrt(15)/5]
    
    # 高斯系数
    GaussA = [5/9, 8/9, 5/9]
    
    # 定义初始区间个数
    n = 1
    
    # 定义初始步长
    h = b - a
    
    # 定义被积函数
    f = lambda x: 0.1*x**2 + 0.2*x - 2*np.sin(x)*np.cos(x)
    
    # 计算准确的积分值
    fun_1, _ = quad(f, a, b)
    
    # 规定精度为1e-20
    error = 1e-20
    
    # 计算初始的积分近似值
    fun_2 = 0
    for i in range(3):
        fun_2 += ((b-a)/2)*GaussA[i]*f((a+b)/2+(b-a)*GaussP[i]/2)
    
    # 当积分近似值与准确值的差的绝对值大于精度时，继续循环
    while abs(fun_1 - fun_2) > error:
        # 增加等分区间
        n += 1
        
        # 重新计算步长
        h = (b - a) / n
        
        # 重置积分近似值
        fun_2 = 0
        
        # 对每个子区间进行积分
        for i in range(n):
            x0 = a + i * h
            x1 = x0 + h
            for j in range(3):
                fun_2 += ((x1-x0)/2) * GaussA[j] * f((x1+x0)/2 + (x1-x0)*GaussP[j]/2)
    
    # 返回积分近似值和等分区间数
    return fun_2, n

# 测试函数
print(gauss_legendre_integr(0, 2*np.pi))
