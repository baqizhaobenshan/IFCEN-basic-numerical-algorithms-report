import numpy as np
from scipy.linalg import norm

def Newton_k(x0, y0, eps=1e-6):
    tol = 1
    r = np.array([x0, y0])
    R = []
    while tol > eps:
        Fx = np.array([r[1] + r[0]**2 - 0.5 - r[0], r[0]**2 - 5*r[0]*r[1] - r[1]])
        R.append([r.tolist(), Fx.tolist()])
        Fx_jacobi = np.array([[2*r[0]-1, 1], [2*r[0]-5*r[1], -5*r[0]-1]])
        rr = r - np.linalg.solve(Fx_jacobi, Fx)
        tol = norm(rr - r)
        r = rr
    return r, Fx, R

# 测试代码
x0, y0 = 1, 0
r, Fx, R = Newton_k(x0, y0)
print("r:\n", r)
print("Fx:\n", Fx)
print("R matrix:")
for i, item in enumerate(R):
    print(f"迭代 {i+1}:")
    print("r: ", item[0])
    print("Fx: ", item[1])