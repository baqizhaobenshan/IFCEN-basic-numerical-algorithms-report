function [ x ] = Chase_method( A, b )
%Chase method 追赶法求三对角矩阵的解
%   A为三对角矩阵的系数，b为等式右端的常数项，返回值x即为最终的解
%   注：A尽量为方阵，b一定要为列向量
%% 求追赶法所需L及U
T = A;
for i = 2 : size(T,1)
    T(i,i-1) = T(i,i-1)/T(i-1,i-1);
    T(i,i) = T(i,i) - T(i-1,i) * T(i,i-1);
end
L = zeros(size(T));
L(logical(eye(size(T)))) = 1;   %对角线赋值1
for i = 2:size(T,1)
    for j = i-1:size(T,1)
        L(i,j) = T(i,j);
        break;
    end
end
U = zeros(size(T));
U(logical(eye(size(T)))) = T(logical(eye(size(T))));
for i = 1:size(T,1)
    for j = i+1:size(T,1)
        U(i,j) = T(i,j);
        break;
    end
end

y = L\b;
x = U\y;
end

