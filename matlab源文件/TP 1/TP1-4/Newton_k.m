%此函数用于求方程组
function [r,Fx,R]=Newton_k(x0,y0,eps) 
if nargin==2 %如果传入变量只有x0时，则默认精度为1e-6；
    eps=1e-6; 
end 

%定义符号变量x和y
syms x y 
F=[y+x^2-0.5-x,x^2-5*x*y-y]; 
%设定一个求雅可比行列式的矩阵
Fjacobi=jacobian(F,[x,y]); 
tol=1; 
r=[x0,y0]; 
R=[]; 
while tol>eps  
    % 通过subs给F函数中的x和y附具体数值
    Fx=subs(F,[x,y],r); 
    Fx=double(Fx); 
    R=[R;[r,Fx]];  %记录x，y以及迭代点的两个函数值
    Fxjacobi=double(subs(Fjacobi,[x,y],r)); %给雅可比行列式赋值
    rr=(r'-Fxjacobi\Fx')';%计算雅可比矩阵Fxjacobi的逆矩阵与Fx的转置矩阵的乘积
    tol=norm(rr-r); 
    r=rr;    
end 
end
