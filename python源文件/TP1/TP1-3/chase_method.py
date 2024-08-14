import numpy as np

def crout(A, B):
    # 获取矩阵A的维度
    n = A.shape[0]
    
    # 初始化L为n*n的零矩阵，U为n*n的单位矩阵
    L = np.zeros((n, n))
    U = np.eye(n)
    
    # 初始化Y和X为n维的零向量
    Y = np.zeros(n)
    X = np.zeros(n)
    
    # 确保A是方阵
    assert A.shape[0] == A.shape[1], "Matrix A must be square"
    
    # 确保A是非奇异的
    assert np.linalg.matrix_rank(A) == n, "Matrix A must be non-singular"
    
    # 对于每个j，更新L和U的值
    for j in range(n):
        L[j:, j] = A[j:, j] - L[j:, :j] @ U[:j, j]
        U[j, j+1:] = (A[j, j+1:] - L[j, :j] @ U[:j, j+1:]) / L[j, j]
    
    # 解线性方程组Ly=b，得到Y
    Y = np.linalg.solve(L, B)
    
    # 解线性方程组Ux=y，得到X
    X = np.linalg.solve(U, Y)
    
    # 返回L, U, Y, X
    return L, U, Y, X
# 定义矩阵A和向量B
A = np.array([
    [10, 3, 0, 0, 0],
    [2, 8, 3, 0, 0],
    [0, 2, 9, 6, 0],
    [0, 0, 7, 11, 9],
    [0, 0, 0, 2, 10]
])
B = np.array([53, 48, 59.4, 86.1, 90.1])

# 调用crout函数
L, U, Y, X = crout(A, B)

# 打印结果
print("L:\n", L)
print("U:\n", U)
print("Y:\n", Y)
print("X:\n", X)