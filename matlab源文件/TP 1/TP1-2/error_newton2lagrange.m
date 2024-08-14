clc,clear;
%在0-2pi内均匀分出41个离散点
x_0=linspace(0,2*pi,41);
y_0=sin(x_0);
%重构0-2pi内均匀的101个离散点
y1=zeros(1,101);
x_a=linspace(0,2*pi,101);
%拉格朗日插值
for i=1:length(x_a)
    y1(i)=lagrange_interp(x_0,y_0,x_a(i));
end
%牛顿插值
y2=newtonzi(x_0,y_0,x_a);
%误差分析
error1=mean((y1-sin(x_a)).^2);
error2=mean((y2-sin(x_a)).^2);
disp("拉格朗日插值误差为"+error1+",牛顿插值误差为"+error2);