## TP 1

### TP 1-1

#### 离散线性方程

**牛顿插值**

若有 $x_0 ,x_1 \dots x_k$以及 $y_0 , y_1 \dots y_k$，则牛顿插值的公式如下
$P_n(x) = f(x_0) + f(x_0,x_1)(x-x_0) + f(x_0,x_1,x_2)(x-x_0)(x-x_1) + \cdots + f(x_0,x_1,\ldots,x_n)(x-x_0)(x-x_1)\cdots(x-x_{n-1})$
其中

$$
\begin{cases}
f(x_0) = f(x_0) \\
f(x_0,x_1) = \frac{f(x_1) - f(x_0)}{x_1 - x_0} \\
f(x_0,x_1,x_2) = \frac{f(x_1,x_2) - f(x_0,x_1)}{x_2 - x_0} \\
\vdots \\
f(x_0,x_1,\ldots,x_k) = \frac{f(x_1,\ldots,x_k) - f(x_0,\ldots,x_{k-1})}{x_k - x_0}
\end{cases}
$$

**拉格朗日插值**

$L_1(x),L_2(x)\dots L_k(x)$ 为拉格朗日插值法的基函数，其公式如下

$$
\begin{cases}
L_1(x) = \frac{(x - x_2)(x - x_3)(x - x_4)...(x - x_n)}{(x_1 - x_2)(x_1 - x_3)(x_1 - x_4)...(x_1 - x_n)} \\
L_2(x) = \frac{(x - x_1)(x - x_3)(x - x_4)...(x - x_n)}{(x_2 - x_1)(x_2 - x_3)(x_2 - x_4)...(x_2 - x_n)} \\
L_3(x) = \frac{(x - x_1)(x - x_2)(x - x_4)...(x - x_n)}{(x_3 - x_1)(x_3 - x_2)(x_3 - x_4)...(x_3 - x_n)} \\
\vdots \\
L_k(x) = \frac{(x - x_1)(x - x_2)(x - x_3)...(x - x_{k-1})(x - x_{k+1})...(x - x_n)}{(x_k - x_1)(x_k - x_2)(x_k - x_3)...(x_k - x_{k-1})(x_k - x_{k+1})...(x_k - x_n)}
\end{cases}
$$

插值的多项式函数为 $P(x) = \sum_{i=1}^n y_i L_i(x)$
</br><br>

#### 流程图

**牛顿插值**

```mermaid
graph LR;
    A[输入已知数据点的x坐标和y坐标] --> B[初始化差商表为离散点的y值];
    B --> C[计算差商表中的值];
    C --> D[输入插值点的x坐标];
    D --> E[计算插值点的函数值];
    E --> F[输出插值点的函数值];
    C --> G[迭代计算差商表中的值];
    G --> C;
```

**拉格朗日插值**

```mermaid
graph LR;
A[输入观测点的x坐标和y坐标] --> B[输入要插值的点的x坐标];
B --> C[初始化插值结果为0];
C --> D[循环计算每个基函数的值];
D --> E[计算当前基函数的值];
E --> F[求和更新插值结果];
F --> G{是否计算完所有基函数};
G -- 是 --> D;
G -- 否 --> H[输出插值结果];
```

</br><br>

#### 代码

**牛顿插值**

```matlab
function y = newtonzi(X, Y, x)
% 牛顿插值函数
% X: 已知数据点的x坐标
% Y: 已知数据点的y坐标
% x: 插值点的x坐标
% y: 插值点的函数值
% 传入离散值的x坐标数组和y坐标数组，以及插值点的x坐标数组，返回插值点的函数值数组

n = length(X); % 已知数据点的数量
A = zeros(n, n); % 初始化差商表
A(:, 1) = Y'; % 第一列为已知数据点的y坐标
for j = 2:n % 遍历差商表的每一列
    for k = j:n % 遍历差商表的每一行
        A(k, j) = (A(k, j-1) - A(k-1, j-1)) / (X(k) - X(k-j+1)); % 计算差商表中的值
    end
end

m = length(x); % 插值点的数量
y = zeros(1, m); % 初始化插值点的函数值
for t = 1:m % 遍历每个插值点
    z = x(t); % 当前插值点的x坐标
    s = Y(1); % 初始化插值点的函数值
    for k = 1:n % 遍历差商表的每一行
        p = 1.0; % 初始化一个乘积
        for j = 1:k-1 % 遍历每个已知数据点
            p = p * (z - X(j)); % 计算乘积
        end
        s = s + A(k, k) * p; % 计算插值点的函数值
    end
    y(t) = s; % 存储插值点的函数值
end
```

**拉格朗日插值**

```matlab
% 拉格朗日插值函数
function y = lagrange(x, y, x0)
    n = length(x);
    y0 = 0;
    for i = 1:n
        L = 1;
        for j = 1:n
            if j ~= i
                L = L .* (x0 - x(j)) / (x(i) - x(j));
            end
        end
        y0 = y0 + y(i) * L;
    end
    y = y0;
end
```

</br><br>

#### 运行结果

**牛顿插值**
![Alt text](1-01.png)

**拉格朗日插值**
![Alt text](1-02.png)

### TP 1-2

#### 代码

```matlab
clc,clear;
%在0-2pi内均匀分出41个离散点
x_0=linspace(0,2*pi,41);
y_0=sin(x_0);
%重构0-2pi内均匀的101个离散点
y1=zeros(1,101);
x_a=linspace(0,2*pi,101);
%拉格朗日插值
for i=1:length(x_a)
    y1(i)=lagrange_interp(x_0,y_0,x_a(i));
end
%牛顿插值
y2=newtonzi(x_0,y_0,x_a);
%误差分析
error1=mean((y1-sin(x_a)).^2);
error2=mean((y2-sin(x_a)).^2);
disp("拉格朗日插值误差为"+error1+",牛顿插值误差为"+error2);
```


#### 运行结果以及结果分析

![Alt text](1-03.png)

### TP 1-3

#### 离散线性方程

追赶法中的 Crout 分解，以下为 Crout 分解矩阵的分解过程

$$
\begin{aligned}
&\begin{bmatrix}L_{11}&0&0\\L_{21}&L_{22}&0\\L_{31}&L_{32}&L_{33}\end{bmatrix}*\begin{bmatrix}1&U_{12}&U_{13}\\0&1&U_{23}\\0&0&1\end{bmatrix}=\begin{bmatrix}L_{11}&L_{11}*U_{12}&L_{11}*U_{13}\\L_{21}&L_{21}*U_{12}+L_{22}&L_{21}*U_{13}+L_{22}*U_{23}\\L_{31}&L_{31}*U_{12}+L_{32}&L_{31}*U_{13}+L_{32}*U_{23}+L_{33}\end{bmatrix} \\
\\
&\begin{bmatrix}L_{11}&L_{11}*U_{12}&L_{11}*U_{13}\\L_{21}&L_{21}*U_{12}+L_{22}&L_{21}*U_{13}+L_{22}*U_{23}\\L_{31}&L_{31}*U_{12}+L_{32}&L_{31}*U_{13}+L_{32}*U_{23}+L_{33}\end{bmatrix}=\begin{bmatrix}A_{11}&A_{12}&A_{13}\\A_{21}&A_{22}&A_{23}\\A_{31}&A_{32}&A_{33}\end{bmatrix}
\end{aligned}
$$

所以有

$$
\begin{aligned}L_{ij}&=A_{ij}-\sum_{k=1}^{j-1}L_{ik}*U_{kj}\\
\\U_{ij}&=\frac{A_{ij}-\sum_{k=1}^{i-1}L_{ik}*U_{kj}}{U_{ii}}\end{aligned}
$$

#### 流程图

```mermaid
graph LR
A[开始] --> B{"迭代j=1:n"}
B --> |是| C["更新L"]
C --> D["更新U"]
D --> B
B --> |否| E["解下三角线性方程组Ly=b"]
E --> F["解上三角线性方程组Ux=y"]
F --> G[结束]
```

#### 代码

```matlab
function [L, U,Y,X] = crout(A,B)
  [m, n] = size(A);
  L = zeros(n, n);%初始化零矩阵
  U = eye(n);%初始化单位矩阵
  Y=zeros(1,n);
  X=zeros(1,n);
  assert(m == n)%断言矩阵为方阵
  assert(length(eigs(A)) == n)%断言矩阵为非奇异的
  for j = 1:n
    L(j:n, j) = A(j:n, j) - L(j:n, 1:j-1) * U(1:j-1, j);%更新矩阵中第j到n行的列
    U(j, j+1:n) = 1 / L(j, j) * (A(j, j+1:n) - L(j, 1:j-1) * U(1:j-1, j+1:n));%更新矩阵中第j+1到n列的行
  end
Y = L\B;%%解Ly=b
X = U\Y;%%解Ux=y
end
```


#### 运行结果以及结果分析


![Alt text](Snipaste_2023-11-20_18-26-26.png)
![Alt text](Snipaste_2023-11-20_18-27-00.png)


### TP 1-4

#### 离散线性方程

$$
\begin{aligned}
&x_{n+1}=x_n-\left[J_F(x_n)\right]^{-1}F(x_n)\\
&\left.F=\left[\begin{array}{c}y+x^2-0.5-x\\x^2-5xy-y\end{array}\right.\right]=\left[\begin{array}{c}-0.5\\1\end{array}\right]=,F’=\left[\begin{array}{cc}2x-1&1\\2x-5y&-5x-1\end{array}\right]=\left[\begin{array}{cc}1&1\\2&-6\end{array}\right] \\
&X_1=\begin{bmatrix}1\\0\end{bmatrix}-\begin{bmatrix}1&1\\2&-6\end{bmatrix}^{-1}\begin{bmatrix}-0.5\\1\end{bmatrix}=\begin{bmatrix}1.25\\0.25\end{bmatrix} \\
&F=\begin{bmatrix}0.0625\\-0.25\end{bmatrix}=,F'=\begin{bmatrix}1.5&1\\1.25&-7.25\end{bmatrix} \\
&\left.X_2=\left[\begin{matrix}1.25\\0.25\end{matrix}\right.\right]-\left[\begin{matrix}1.5&1\\1.25&-7.25\end{matrix}\right]^{-1}\left[\begin{matrix}0.0625\\-0.25\end{matrix}\right]=\left[\begin{matrix}1.2332\\0.2126\end{matrix}\right] \\
&...
\end{aligned}
$$

#### 流程图

```mermaid
graph LR;
    A[开始] --> B[输入初始值x0];
    B --> C[定义符号变量x和y];
    C --> D["计算方程组F(x,y)"];
    D --> E[计算雅可比矩阵Fjacobi];
    E --> F[计算新的迭代点r];
    F --> G[判断是否满足收敛条件];
    G -- 是 --> H[输出迭代结果];
    G -- 否 --> I[更新迭代点的值r];
    I --> D;
```

#### 代码

```matlab
%此函数用于求方程组
function [r,Fx,R]=Newton_k(x0,y0,eps)%定义两个函数
if nargin==2 %如果传入变量只有x0时，则默认精度为1e-6；
    eps=1e-6;
end

%定义符号变量x和y
syms x y
F=[y+x^2-0.5-x,x^2-5*x*y-y];
%设定一个求雅可比行列式的矩阵
Fjacobi=jacobian(F,[x,y]);
tol=1;
r=[x0,y0];
R=[];
while tol>eps
    % 通过subs给F函数中的x和y附具体数值
    Fx=subs(F,[x,y],r);
    Fx=double(Fx);
    R=[R;[r,Fx]];  %记录x，y以及迭代点的两个函数值
    Fxjacobi=double(subs(Fjacobi,[x,y],r)); %给雅可比行列式赋值
    rr=(r'-Fxjacobi\Fx')';%计算雅可比矩阵Fxjacobi的逆矩阵与Fx的转置矩阵的乘积
    tol=norm(rr-r);
    r=rr;
end
end
```

#### 运行结果以及结果分析
![Alt text](Snipaste_2023-11-20_22-53-15.png)

### TP 2-1

#### 离散线性方程

**单点积分格式**

$$
\int_a^b f(x)dx \approx (b-a) \cdot f(a)
$$

**梯形积分格式**

$$
\int_a^b f(x)dx \approx \frac{b-a}{2} \cdot (f(a)+f(b))
$$

**高斯积分格式**
在区间$[-1,1]$上，高斯-勒让德求积公式为

$$
\int_{-1}^{1} f(x)dx \approx \sum_{k=0}^{n} A_kf(x_k)
$$

三点高斯-勒让德求积公式为

$$
\int_{-1}^{1} f(x) dx \approx \frac{5}{9}f(\frac{-\sqrt{15}}{5})+\frac{8}{9}f(0)+\frac{5}{9}f(\frac{-\sqrt{15}}{5})
$$

若积分区间不是$[-1,1]$而是$[a,b]$时
则

$$
\int_{a}^{b} f(x) dx=\frac{b-a}{2}\int_{-1}^{1}f(\frac{b-a}{2}x+\frac{a-b}{2})dx=\frac{b-a}{2}\sum_{N}^{k=1}\omega_{N,k}f(\frac{a+b}{2}+\frac{b-a}{2}x_{N,k})
$$

#### 流程图以及实现思路

对于三种积分方式，虽然计算公式不同，但是实现思路是一样的

```mermaid
graph LR
A(开始)-->B(设置积分函数)
B-->C("初始化等分区间为1，并计算准确积分")
C-->D(初始化步长)
D-->F(计算近似的积分值)
F-->G(比较准确积分与数值积分以判断误差是否满足要求)
G--是-->H(输出近似的积分值和等分区间数)
G--否-->D
H-->I(结束)
```

#### 代码

**单点积分**

```matlab
function [intgr,num] = single_point_integr(a,b)%a,b为积分区间
f=@(x) 0.1.*x.^2+0.2.*x-2.*sin(x).*cos(x);%设置积分函数
n=1;%初始等分区间为1
h=b-a;%初始的步长
tol=1e-6;%设置精度
fun_1=integral(f,a,b);%matlab计算出的准确的积分值
fun_2=f(b)*h; %调用方程函数
while abs(fun_2-fun_1)>tol
 n=n+1;
 h=(b-a)/n;
 fun_2=0;
    for i=0:n-1
          x=a+i*h;
          x1=x+h;
          fun_2=fun_2+h*f(x1);%按区间累加以获得[a,b]区间的数值积分
    end
    if n>9999%若划分的区间很大时仍无法达到所需精度则跳出循环
        break;
    end
end
intgr=fun_2;num=n;
end
```

**梯形积分**

```matlab
function [intgr,num] = td_integr(a,b)%a,b为积分区间
f=@(x) 0.1.*x.^2+0.2.*x-2.*sin(x).*cos(x);%设置积分函数
n=1;%初始等分区间为1
h=b-a;%初始的步长
tol=1e-6;
fun_1=integral(f,a,b);%matlab计算出的准确的积分值
fun_2=(f(a)+f(b))*(h/2); %调用方程函数
while abs(fun_2-fun_1)>tol
 n=n+1;
 h=(b-a)/n;
 fun_2=0;
    for i=0:n-1
          x=a+i*h;
          x1=x+h;
          fun_2=fun_2+(h/2)*(f(x)+f(x1));%按区间累加以获得[a,b]区间的数值积分
    end
end
intgr=fun_2;num=n;
end
```

**三点高斯-勒让德积分**

```matlab
function [result,num]= gauss_legendre_integr(a,b)
    % 高斯-勒让德求积公式三点情况
    GaussP = [-sqrt(15)/5, 0, sqrt(15)/5]; % 高斯点
    GaussA = [5/9,8/9,5/9]; % 高斯系数
    n=1;%定义初始区间个数
    h=b-a;%定义初始步长
    f=@(x) 0.1.*x.^2+0.2.*x-2.*sin(x).*cos(x);
    fun_1 = integral(f, a, b); % 计算准确的积分值
    error=1e-20;%%规定精度为1e-20
    fun_2 =0;
    for i=1:3
        fun_2=((b-a)/2)*GaussA(i)*f((a+b)/2+(b-a)*GaussP(i)/2)+fun_2;
    end
    while abs(fun_1-fun_2) > error
     n=n+1;
     h=(b-a)/n;
     fun_2=0;
        for i=0:n-1
          x0=a+i*h;
          x1=x0+h;
          fun_2=fun_2+((x1-x0)/2)*(GaussA(1)*f((x1+x0)/2+(x1-x0)*GaussP(1)/2)+GaussA(2)*f((x1+x0)/2+(x1-x0)*GaussP(2)/2)+GaussA(3)*f((x1+x0)/2+(x1-x0)*GaussP(3)/2));
        end
    end
    result=fun_2;
    num=n;
end

```

#### 运行结果以及结果分析

**单点积分**
![Alt text](Snipaste_2023-11-22_13-00-56.png)
未达到 1e-6 精度的误差，划分区间超过 10000 后退出循环

**梯形积分**
![Alt text](Snipaste_2023-11-22_00-31-16.png)

计算所得的数值积分值为 12.2162，划分出 2034 个区间

**三点高斯-勒让德积分**
![Alt text](Snipaste_2023-11-22_15-52-26.png)
可见达到 1e-20 精度时只划分了 7 个区间

```mermaid
graph LR
A(开始)-->B("设置参数，初始化n、C数组")
B-->C("列出方程组")
C-->E("判断i是否达到数组的数目")
E--是-->D("数值求解n(i)、C(i)(选取不同的求解方法)")
D-->F(绘图)
E--否-->F(绘图)
F-->I(结束)
```

```mermaid
graph LR
    A[初始化] --> B[初始化矩阵T，确定各个参数]
    B --> C[使用二阶差分法更新矩阵T]
    C --> G["比较新矩阵与旧矩阵，T_new-T_old的精度是否达到标准"]
    G --是--> H[跳出循环]
    G--否-->C
    H --> J[将T绘制为热图]
```
### TP 2-2

#### 离散线性方程

**显式法**：

$$
y_{i+1}=y_i+hf(x_i,y_i)
$$

**隐式法**：

$$
y_{i+1}=y_i+hf(x_i+1,y_i+1)
$$

**四阶 R-K 法**：

$$
\begin{cases}
K_1=h \cdot f(t_n,y_n) \\
K_2=h \cdot f(t_n+\frac{1}{2}h,y_n+\frac{1}{2}hK_1) \\
K_3=h \cdot f(t_n+\frac{1}{2}h,y_n+\frac{1}{2}hK_2)\\
K_4=h \cdot f(t_n+h,y_n+hK_3)\\
y_{n+1}=y_n+\frac{1}{6}K_1+\frac{1}{3}K_2+\frac{1}{3}K_3+\frac{1}{6}K_4
\end{cases}
$$

#### 流程图

```mermaid
graph LR
A(开始)-->B("设置参数，初始化n、C数组")
B-->C("列出方程组")
C-->E("判断i是否达到数组的数目")
E--是-->D("数值求解n(i)、C(i)(选取不同的求解方法)")
D-->F(绘图)
E--否-->F(绘图)
F-->I(结束)
```



#### 代码实现

**显式法**

```matlab
function [n,C]=xianshi(h)
rho=0.0022;Lambda=10^(-3);beta=0.0065;lambda=0.078;%定义参数
t=[0 : h :1];%定义时间范围
n=zeros(size(t));
C=zeros(size(t));
len = length(t);
function n1 =Fun1(t,n,C)%定义方程组
n1=(rho-beta).*n/Lambda+lambda.*C;
end
function C1=Fun2(t,n,C)
C1=beta.*n./Lambda-lambda.*C;
end
%定义初值
n(1)=1;      
C(1)=1;      
for i=2:len 
n(i)=n(i-1)+h*Fun1(t(i-1),n(i-1),C(i-1));% 更新下一个dn/dt
C(i)=C(i-1)+h*Fun2(t(i-1),n(i-1),C(i-1));% 更新下一个dC/dt
end
hold on;
grid on;
plot(t, n,'rh',t,C,'bo',LineWidth=1,MarkerSize=3);		
xlabel('t/s');
legend('$n(t)$','$C(t)$');
title("显式法求解中子动力学方程",Interpreter="none");
fig=gcf;
fig.Children(1).Interpreter='latex';
fig.Children(1).Title.Interpreter='latex';
fig.Children(1).FontSize=13.5;
fig.Children(2).FontSize=15;
end
```

<div style="page-break-after: always;"></div>

**隐式法**

```matlab
function [n,C]=yinshi(h)
rho=0.0022;Lambda=10^(-3);beta=0.0065;lambda=0.078;%定义参数
t=[0 : h :1];%定义时间范围
n=zeros(size(t));
C=zeros(size(t));
len = length(t);
function n1 =Fun1(t,n,C)%定义方程组
n1=(rho-beta).*n/Lambda+lambda.*C;
end
function C1=Fun2(t,n,C)
C1=beta.*n./Lambda-lambda.*C;
end
%定义初值
n(1)=1;      
C(1)=1;      
for i=2:len 
A = [1-h*(rho-beta)/Lambda, -h*lambda; -h*beta/Lambda, 1+h*lambda];
B = [n(i-1); C(i-1)];
%求解线性方程组
X = linsolve(A,B);
% 提取n(i)和C(i)
n(i) = X(1);
C(i) = X(2);
end
% hold on;
% grid on;
% plot(t, n,'rh',t,C,'bo',LineWidth=1,MarkerSize=3);		
% xlabel('t/s');
% legend('$n(t)$','$C(t)$');
% title("隐式法求解中子动力学方程",Interpreter="none");
% fig=gcf;
% fig.Children(1).Interpreter='latex';
% fig.Children(1).Title.Interpreter='latex';
% fig.Children(1).FontSize=13.5;
% fig.Children(2).FontSize=15;
end
```

<div style="page-break-after: always;"></div>

**四阶 R-K 法**

```matlab
function [n,C]=Four_RK(h)
rho=0.0022;Lambda=10^(-3);beta=0.0065;lambda=0.078;%定义参数
t=[0 : h :1];%定义时间范围
n=zeros(size(t));
C=zeros(size(t));
len = length(t);
function n1 =Fun1(t,n,C)%定义方程组
n1=(rho-beta).*n/Lambda+lambda.*C;
end
function C1=Fun2(t,n,C)
C1=beta.*n./Lambda-lambda.*C;
end
%定义初值
n(1)=1;      
C(1)=1;      
for i=2:len
%Kmn 表示第n次迭代，以m次导数作为f(x,y)
K11 = Fun1(t(i-1),n(i-1),C(i-1));  
K21 = Fun2(t(i-1),n(i-1),C(i-1));
K12 = Fun1(t(i-1)+1/2*h , n(i-1)+1/2*h*K11 , C(i-1)+1/2*h*K21);
K22 = Fun2(t(i-1)+1/2*h , n(i-1)+1/2*h*K11 , C(i-1)+1/2*h*K21);
K13 = Fun1(t(i-1)+1/2*h , n(i-1)+1/2*h*K12 , C(i-1)+1/2*h*K22);
K23 = Fun2(t(i-1)+1/2*h , n(i-1)+1/2*h*K12 , C(i-1)+1/2*h*K22);
K14 = Fun1(t(i-1)+h , n(i-1)+h*K13 , C(i-1)+h*K23);
K24 = Fun2(t(i-1)+h , n(i-1)+h*K13 , C(i-1)+h*K23);
n(i) = n(i-1)+h/6*(K11 + 2*K12 + 2*K13 + K14);           % 更新下一个dn/dt
C(i) = C(i-1)+h/6*(K21 + 2*K22 + 2*K23 + K24);         % 更新下一个dC/dt
end
hold on;
grid on;
plot(t, n,'rh',t,C,'bo',LineWidth=1,MarkerSize=3);		
xlabel('t/s');
legend('$n(t)$','$C(t)$');
title("四阶R-K法求解中子动力学方程",Interpreter="none");
fig=gcf;
fig.Children(1).Interpreter='latex';
fig.Children(1).Title.Interpreter='latex';
fig.Children(1).FontSize=13.5;
fig.Children(2).FontSize=15;
end
```



#### 结果分析

测试代码：`result_test.m`

```matlab
clc,clear;
h_values = [0.1, 0.01, 0.001,0.0001,0.00001];  % 定义步长数组

for i = 1:length(h_values)
    h = h_values(i);
    
    [n, C] = Four_RK(h);
    fprintf('Four_RK - Step size: %f, Last n: %f, Last C: %f\n', h, n(end), C(end));
    
    [n, C] = xianshi(h);
    fprintf('xianshi - Step size: %f, Last n: %f, Last C: %f\n', h, n(end), C(end));
    
    [n, C] = yinshi(h);
    fprintf('yinshi - Step size: %f, Last n: %f, Last C: %f\n', h, n(end), C(end));
    fprintf("------------------------------------------------\n");
end
```
测试代码: `stability_test.m`
```matlab 
clc,clear;
% 定义一个步长数组
h = [0.1, 0.05, 0.01, 0.005, 0.001];
% 定义一个结果数组
result = zeros(2, length(h));
% 定义一个时间数组
time = zeros(1, length(h));
% 定义一个参考解，用最小的步长求解
[n_ref, C_ref] = Four_RK(h(end));
% 循环遍历不同的步长
for i = 1:length(h)
    % 记录开始时间
    tic;
    % 调用隐式、显式或四阶R-K函数求解
    [n, C] =Four_RK(h(i));
    % 记录结束时间
    toc;
    % 计算运行时间
    time(i) = toc - tic;
    % 计算结果的差异
    result(1, i) = abs(n(end) - n_ref(end));
    result(2, i) = abs(C(end) - C_ref(end));
end
% 显示结果
% 显示结果
disp('不同步长的n的差异：');
disp(result(1,:));
disp('不同步长的C的差异：');
disp(result(2,:));
```

首先比较三种不同算法的结果，可以看见四阶 R-K 法收敛得较快：
![alt text](微分方程结果对比.png)

然后是稳定性与速度对比，比较三种不同的求解方法的原理为记录每种方法下不同步长计算所需的时长，并且定义步长最小时计算出来的结果为稳定解，比较不同步长下计算出来的结果，从而评估算法的稳定性以及速度.

隐式:
![alt text](yinshi.png)
显式:
![alt text](xianshi.png)
四阶R-K法:
![alt text](4-RK.png)

1. 我们可以发现三者中显式法的计算耗时最少，四阶 R-K 法的耗时要比隐式法的耗时短，这可能是因为在隐式法求解时，使用了 matlab 内置的 `linsolve` 函数来线性方程组，这导致了当步长较小时，由于循环的原因，隐式法的耗时要比四阶 R-K 法要更长；

2. 四阶 R-K 法中 n、C 差异是最小的，说明四阶 R-K 法的稳定性也是最好的，收敛速度较快；

![alt text](四阶R-K法.jpg)
## TP 3

### TP 3-1

#### 离散线性方程

采用二阶中心差分格式，即：\\

$$
\frac{\partial^2 T}{\partial x^2}=\frac{T_{i+1}+T_{i-1}-2Ti}{(\Delta x)^2}
$$

$$
\frac{\partial^2 T}{\partial y^2}=\frac{T_{i+1}+T_{i-1}-2Ti}{(\Delta y)^2}
$$

可以推导出

$$
T(i,j)=\frac{\frac{T(i+1,j)+T(i-1,j)}{(dx)^2}+\frac{T(i,j+1)+T(i,j-1)}{(dy)^2}-4.5\cdot e^{1.5(dx\cdot i+dy\cdot j)}}{\frac{2}{dx^2}+\frac{2}{dy^2}}
$$

#### 流程图

```mermaid
graph LR
    A[初始化] --> B[初始化矩阵T，确定各个参数]
    B --> C[使用二阶差分法更新矩阵T]
    C --> G["比较新矩阵与旧矩阵，T_new-T_old的精度是否达到标准"]
    G --是--> H[跳出循环]
    G--否-->C
    H --> J[将T绘制为热图]
```



#### 代码实现

```matlab
clc,clear;
Nx=100;Ny=100;%我们对x轴分成Nx段，有Nx+1个点，对t轴分成Nt段，有Nt+1个点
xd=1;yd=1;
dx=xd/Nx;dy=yd/Ny;
T=zeros(Nx+1,Ny+1);%生成(Nx+1)*(Ny+1)的矩阵
T(1,:) = 293 + exp(3/2*linspace(0,xd,Nx+1));%T(0,y)=293+exp(3y/2)
T(end,:)=293+exp((3/2)*(1+linspace(0,xd,Nx+1)));
T(:,1)=293 + exp(3/2*linspace(0,yd,Ny+1));
T(:,end)=293 + exp((3/2)*(1+linspace(0,yd,Ny+1)));
error = 1e-3; % 误差阈值
max_iterations = 100000; % 最大迭代次数
% 进行迭代计算
T_news=zeros(Nx+1,Ny+1);
for iteration = 1:max_iterations
    % 备份旧的温度矩阵
    T_old = T;
    for i=2:Nx
        for j=2:Ny
            T(i,j)=((T(i+1,j)+T(i-1,j))/(dx)^2+(T(i,j+1)+T(i,j-1))/(dy)^2-4.5*exp(1.5*(dx*i+dy*j)))/(2/dx^2+2/dy^2);
        end
    end
    T_new=T;
    diff = abs(T - T_old);
    if max(diff(:)) < error
        break;
    end
    % 更新温度矩阵
    T = T_new;
end
x = dx*(0:Nx);  
y = dy*(0:Ny);
[X,Y]=meshgrid(x,y);
s=surf(X,Y,T);
s.EdgeColor="none";
zlabel("T/K");
colorbar;
title("有限差分法求解二维传热方程");
fig=gcf;
fig.Children(1).FontSize=13.5;
fig.Children(2).FontSize=15;

T_exact = @(x,y) 293+exp(1.5*(x+y)); 
T_exact_array=T_exact(X,Y);
T_error=T_exact_array-T;
figure;
s=surf(X,Y,T_error);
s.EdgeColor="none";
zlabel("$\Delta T/K$",Interpreter="latex");
colorbar;
title("有限差分法求解二维传热方程的误差");
fig=gcf;
fig.Children(1).FontSize=13.5;
fig.Children(2).FontSize=15;
```



#### 结果分析

使用解析解 $T=293+exp\left(\frac{3}{2}(x+y)\right)$ 计算得到准确值 $T_{exact}$，$\Delta T = T_{exact} - T$，数值解曲面图 T 在左，误差 $\Delta T$ 的曲面图在右。

| 精细度与迭代次数 | 数值解曲面图 | 误差曲面图 |
| ------------ | ------------- | ------------ |
| $11 \times 11$，迭代次数 109 次 | ![数值解曲面图](导热方程步长11_11.png) | ![误差曲面图](导热方程步长11_11误差.png) |
| $51 \times 51$，迭代次数 1914 次 | ![数值解曲面图](导热方程步长51_51.png) | ![误差曲面图](导热方程步长51_51误差.png) |
| $101 \times 101$，迭代次数 6249 次 | ![数值解曲面图](导热方程步长101_101.png) | ![误差曲面图](导热方程步长101_101误差.png) |

可以发现：
1. 解析解与数值解的主要误差来源于网格中间的区域；
2. 在三者中，网格精细程度居中的 $51 \times 51$ 网格误差最小，说明网格越精细，结果不一定越准确。
