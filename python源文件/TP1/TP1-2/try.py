import numpy as np
from ..TP1-1.newtonzi import newton_interp
from ..TP1-1.lagrange import lagrange

# 在0 - 2pi内均匀分出41个离散点
x_0 = np.linspace(0, 2 * np.pi, 41)
y_0 = np.sin(x_0)

# 重构0 - 2pi内均匀的101个离散点
y1 = np.zeros(101)
x_a = np.linspace(0, 2 * np.pi, 101)

# 拉格朗日插值
for i in range(len(x_a)):
    y1[i] = lagrange(x_0, y_0, x_a[i])

# 牛顿插值
y2 = newton_interp(x_0, y_0, x_a)

# 误差分析
error1 = np.mean((y1 - np.sin(x_a)) ** 2)
error2 = np.mean((y2 - np.sin(x_a)) ** 2)

print(f"拉格朗日插值误差为 {error1}, 牛顿插值误差为 {error2}")
