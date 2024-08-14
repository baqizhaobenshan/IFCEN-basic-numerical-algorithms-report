function y = newtonzi(X, Y, x)
% 牛顿插值函数
% X: 已知数据点的x坐标
% Y: 已知数据点的y坐标
% x: 插值点的x坐标
% y: 插值点的函数值

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