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