import time
import four_R_K
import xianshi
import yinshi

# 定义参数
rho = 0.0022
Lambda = 10**(-3)
beta = 0.0065
lambda_ = 0.078

# 定义步长列表
h_values = [0.1, 0.01, 0.001, 0.0001]

# 定义求解方法列表
methods = [xianshi.solve_explicit, yinshi.solve_implicit, four_R_K.solve_RK4]  # 这里假设你已经定义了这三个函数

# 对每种方法进行测试
for method in methods:
    print(f"Testing {method.__name__}...")
    for h in h_values:
        start_time = time.time()
        t, n, C = method(rho, Lambda, beta, lambda_, h)
        end_time = time.time()
        print(f"Step size: {h}, Time elapsed: {end_time - start_time} seconds")
        print(f"Final values: n={n[-1]}, C={C[-1]}")
        print()  # 换步长时加入空行
    print('-' * 50)  # 不同方法之间加入分割线