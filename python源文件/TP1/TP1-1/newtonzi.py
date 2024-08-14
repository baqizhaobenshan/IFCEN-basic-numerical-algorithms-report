import numpy as np
def newton_diff(X, Y):
    n = len(Y)
    coef = np.zeros([n, n]) # 创建一个n*n的零矩阵
    coef[:,0] = Y # 第一列是Y

    # 计算差商
    for j in range(1,n):
        for i in range(n-j):
            # 差商的计算公式
            coef[i][j] = (coef[i+1][j-1] - coef[i][j-1]) / (X[i+j]-X[i])

    return coef[0, :] # 返回第一行，这是差商的系数

def newton_interp(X, Y, x):
    a = newton_diff(X, Y) # 计算差商的系数
    n = len(a) - 1
    p = a[n]

    # 构建插值多项式
    for k in range(1, n + 1):
        p = a[n - k] + (x - X[n - k]) * p

    return p # 返回插值结果

# 定义已知的数据点
X = np.array([1, 2, 3, 4, 5])
Y = np.array([1, 16, 81, 256,625]) 

x = 3.5 # 需要插值的点

# 计算插值结果
y = newton_interp(X, Y, x)

# 打印插值结果
print(y)
