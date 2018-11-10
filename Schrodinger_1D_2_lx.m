%Schrodinger_1D_2
%2018/11/1 林祥
%-1/2*phin(s)''+v(s)phin(s)=en*phin(s)
% i). v(s)=0,0<=s<=1; v(s)=inf.s<0 or s>1;
% ii). v(s)=-(1-s^2/2);   

%ii). 本征值理论为n-1/2;n=0,1,2,3...n为节点数
clc;  clear;  format long;
%%t=s; y=phi;
h=0.01;   %时间步长
y0=0;  delta=1;    parity=-1;    %奇宇称
%y0=0.3;  delta=0;     parity=1;  %偶宇称
en=5;  den=0.0002;  %本征值的初始打靶值和初始步长
ETol=1E-5;
%----------------求解--------------
[n,tn,t,yn]=Numerov2(h,y0,delta,en);  %振荡区
while abs(yn(n+1))>ETol  %指数衰减区验证
    tn=tn+h; n=n+1; t=[t,tn];
    yn(n+1)=2*(1-5*h^2/12*(2*en+2-(t(n))^2))*yn(n)-(1+h^2/12*(2*en+2-(t(n-1))^2))*yn(n-1);
    yn(n+1)=yn(n+1)/(1+h^2/12*(2*en+2-(t(n+1))^2));
    if abs(yn(n+1))>abs(yn(n))   %若不衰减，改变en,重新计算
        en=en+den;     [n,tn,t,yn]=Numerov2(h,y0,delta,en);  
    end
end
%---------------归一化---------------
sum2=0;  
for i=1:n          %梯形积分公式
    sum2=sum2+h*((yn(i)+yn(i+1))/2)^2;  
end
sum2=2*sum2;   %加上负轴
yn=yn./sqrt(sum2); 
%----------------画图-----------------
figure(1);  set(gca,'Fontsize',16);
plot(t,yn,'r.',-1.*t,parity.*yn,'r.');   hold on; %补全负轴
xlabel('s');  ylabel('\phi');
legend(sprintf('本征值为%.4f的波函数',en));
title('ii). \phi(s)的函数曲线');
grid on;
%---------函数定义----------
function [n,tn,t,y]=Numerov2(h,y0,delta,en)%震荡区波函数求解
    y(1)=y0; y(2)=y0+h*delta; 
    tn=sqrt(2*en+2);    n=floor(tn/h);  t=0:h:tn;%振荡区
    for i=2:n
        y(i+1)=2*(1-5*h^2/12*(2*en+2-(t(i))^2))*y(i)-(1+h^2/12*(2*en+2-(t(i-1))^2))*y(i-1);
        y(i+1)=y(i+1)/(1+h^2/12*(2*en+2-(t(i+1))^2));
    end
end