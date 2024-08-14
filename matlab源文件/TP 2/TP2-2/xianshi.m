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


