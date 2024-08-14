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



