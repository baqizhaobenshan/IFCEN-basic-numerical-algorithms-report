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

