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
