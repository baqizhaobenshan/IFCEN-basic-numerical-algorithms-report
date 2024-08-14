import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
# 定义参数
Nx = 50
Ny = 50
xd = 1
yd = 1
dx = xd / Nx
dy = yd / Ny
T = np.zeros((Nx+1, Ny+1))
T[0, :] = 293 + np.exp(3/2 * np.linspace(0, xd, Nx+1))
T[-1, :] = 293 + np.exp((3/2) * (1 + np.linspace(0, xd, Nx+1)))
T[:, 0] = 293 + np.exp(3/2 * np.linspace(0, yd, Ny+1))
T[:, -1] = 293 + np.exp((3/2) * (1 + np.linspace(0, yd, Ny+1)))
error = 1e-3
max_iterations = 100000

# 进行迭代计算
for iteration in range(max_iterations):
    T_old = T.copy()
    for i in range(1, Nx):
        for j in range(1, Ny):
            T[i, j] = ((T[i+1, j] + T[i-1, j]) / dx**2 +
                       (T[i, j+1] + T[i, j-1]) / dy**2 -
                       4.5 * np.exp(1.5 * (dx*i + dy*j))) / (2/dx**2 + 2/dy**2)
    diff = np.abs(T - T_old)
    if np.max(diff) < error:
        break

# 绘制结果
x = dx * np.arange(Nx+1)
y = dy * np.arange(Ny+1)
X, Y = np.meshgrid(x, y)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 创建一个颜色映射对象
cmap = cm.coolwarm

# 使用颜色映射对象和数据 T 来生成颜色
norm = plt.Normalize(T.min(), T.max())
colors = cmap(norm(T))

# 在创建 surf 对象时设置 facecolors 参数，并使用 rcount 和 ccount 参数来限制生成的多边形的数量
rcount, ccount = T.shape
surf = ax.plot_surface(X, Y, T, facecolors=colors, rcount=rcount, ccount=ccount, linewidth=0.5)

# 创建 colorbar
m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array(T)
fig.colorbar(m)

ax.set_zlabel('T/K')
plt.title('有限差分法求解二维传热方程')
plt.show()
# 计算解析解
T_analytical = 293 + np.exp(1.5 * (X + Y))

# 计算差异
diff = T_analytical - T

# 绘制差异的3D图
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 创建一个颜色映射对象
cmap = cm.coolwarm

# 使用颜色映射对象和数据 diff 来生成颜色
norm = plt.Normalize(diff.min(), diff.max())
colors = cmap(norm(diff))

# 在创建 surf 对象时设置 facecolors 参数，并使用 rcount 和 ccount 参数来限制生成的多边形的数量
rcount, ccount = diff.shape
surf = ax.plot_surface(X, Y, diff, facecolors=colors, rcount=rcount, ccount=ccount, linewidth=0.5)

# 创建 colorbar
m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array(diff)
fig.colorbar(m)

ax.set_zlabel('Difference')
plt.title('差异曲面图')
plt.show()
print("迭代次数：", iteration)