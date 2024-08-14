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
