def lagrange(x, y, x0):
    # 获取已知点的数量
    n = len(x)
    # 初始化插值结果
    y0 = 0
    # 遍历每个已知点
    for i in range(n):
        # 初始化拉格朗日基函数
        L = 1
        # 计算拉格朗日基函数
        for j in range(n):
            if j != i:
                L *= (x0 - x[j]) / (x[i] - x[j])
        # 更新插值结果
        y0 += y[i] * L
    # 返回插值结果
    return y0

# 已知的点
x = [1, 2, 3, 4,5]
y = [1, 16, 81, 256,625]

# 想要插值的点
x0 = 3.5

# 使用拉格朗日插值函数
result = lagrange(x, y, x0)

# 打印插值结果
print(f"The interpolated value at {x0} is {result}")
