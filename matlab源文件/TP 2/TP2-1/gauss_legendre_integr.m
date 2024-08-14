function [result,num]= gauss_legendre_integr(a,b)
    % 高斯-勒让德求积公式三点情况
    GaussP = [-sqrt(15)/5, 0, sqrt(15)/5]; % 高斯点
    GaussA = [5/9,8/9,5/9]; % 高斯系数
    n=1;%定义初始区间个数
    h=b-a;%定义初始步长
    f=@(x) 0.1.*x.^2+0.2.*x-2.*sin(x).*cos(x);
    fun_1 = integral(f, a, b); % 计算准确的积分值
    error=1e-20;
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
