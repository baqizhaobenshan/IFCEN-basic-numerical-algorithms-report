from scipy.integrate import quad
import numpy as np
def td_integr(a, b):
    # 定义被积函数
    f = lambda x: 0.1*x**2 + 0.2*x - 2*np.sin(x)*np.cos(x)

    # 初始化等分区间为1
    n = 1

    # 计算初始步长
    h = b - a

    # 设置精度
    tol = 1e-6

    # 使用 scipy 的 quad 函数计算准确的积分值
    fun_1, _ = quad(f, a, b)

    # 计算初始的积分近似值，使用梯形公式
    fun_2 = (f(a) + f(b)) * (h / 2)

    # 当积分近似值与准确值的差的绝对值大于精度时，继续循环
    while abs(fun_2 - fun_1) > tol:
        # 增加等分区间
        n += 1

        # 重新计算步长
        h = (b - a) / n

        # 重置积分近似值
        fun_2 = 0

        # 对每个子区间进行积分，使用梯形公式
        for i in range(n):
            x = a + i * h
            x1 = x + h
            fun_2 += (h / 2) * (f(x) + f(x1))

        # 如果等分区间过大，跳出循环
        if n > 9999:
            break

    # 返回积分近似值和等分区间数
    return fun_2, n
# 测试函数
print(td_integr(0, 2*np.pi))
