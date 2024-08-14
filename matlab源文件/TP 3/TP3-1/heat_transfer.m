clc,clear;
Nx=100;Ny=100;%我们对x轴分成Nx段，有Nx+1个点，对t轴分成Nt段，有Nt+1个点
xd=1;yd=1;
dx=xd/Nx;dy=yd/Ny;
T=zeros(Nx+1,Ny+1);%生成(Nx+1)*(Ny+1)的矩阵
T(1,:) = 293 + exp(3/2*linspace(0,xd,Nx+1));%T(0,y)=293+exp(3y/2)
T(end,:)=293+exp((3/2)*(1+linspace(0,xd,Nx+1)));
T(:,1)=293 + exp(3/2*linspace(0,yd,Ny+1));
T(:,end)=293 + exp((3/2)*(1+linspace(0,yd,Ny+1)));
error = 1e-3; % 误差阈值
max_iterations = 100000; % 最大迭代次数
% 进行迭代计算
T_news=zeros(Nx+1,Ny+1);
for iteration = 1:max_iterations
    % 备份旧的温度矩阵
    T_old = T;
    for i=2:Nx
        for j=2:Ny
            T(i,j)=((T(i+1,j)+T(i-1,j))/(dx)^2+(T(i,j+1)+T(i,j-1))/(dy)^2-4.5*exp(1.5*(dx*i+dy*j)))/(2/dx^2+2/dy^2);
        end
    end
    T_new=T;
    diff = abs(T - T_old);
    if max(diff(:)) < error
        break;
    end
    % 更新温度矩阵
    T = T_new;
end
x = dx*(0:Nx);  
y = dy*(0:Ny);
[X,Y]=meshgrid(x,y);
s=surf(X,Y,T);
s.EdgeColor="none";
zlabel("T/K");
colorbar;
title("有限差分法求解二维传热方程");
fig=gcf;
fig.Children(1).FontSize=13.5;
fig.Children(2).FontSize=15;

T_exact = @(x,y) 293+exp(1.5*(x+y)); 
T_exact_array=T_exact(X,Y);
T_error=T_exact_array-T;
figure;
s=surf(X,Y,T_error);
s.EdgeColor="none";
zlabel("$\Delta T/K$",Interpreter="latex");
colorbar;
title("有限差分法求解二维传热方程的误差");
fig=gcf;
fig.Children(1).FontSize=13.5;
fig.Children(2).FontSize=15;

