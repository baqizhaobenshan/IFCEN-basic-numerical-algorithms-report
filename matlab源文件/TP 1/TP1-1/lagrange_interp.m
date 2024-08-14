% 拉格朗日插值函数
function y = lagrange_interp(x, y, x0)
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

