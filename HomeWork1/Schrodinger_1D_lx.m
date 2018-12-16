%Schrodinger_1D
%2018/11/1 林祥
%-1/2*phin(s)''+v(s)phin(s)=en*phin(s)
% i). v(s)=0,0<=s<=1; v(s)=inf.s<0 or s>1;
% ii). v(s)=-(1-s^2/2);
clc;  clear;  format long;
%%t=s; y=phi;
h=0.01;  n=1/h;   t=0:h:1;  %时间范围及 步数
y0=0; yn=0;  %边值
delta=1;  %归一化方程条件!?, delta可为任意值
en=20;  den=0.01;  %本征值的初始打靶值和初始步长
ETol=1E-5;
%----------------求解--------------
yn1=Numerov(n,h,y0,delta,en);  %tspan=[0,1]
while   abs(den)>ETol
    en=en+den;
    yn2=Numerov(n,h,y0,delta,en);
    if   (yn2(n+1))*(yn1(n+1))>0  %对en的搜索实现yn=0的边值条件
    else
        en=en-den;
        den=den/2;
    end
end
%---------------归一化---------------
sum1=0;  
for i=1:n          %梯形积分公式
    sum1=sum1+h*((yn2(i)+yn2(i+1))/2)^2;  
end
yn2=yn2./sqrt(sum1);    
%----------------画图-----------------
figure(1);  set(gca,'Fontsize',16);
plot(t,yn2,'r.');   hold on;
xlabel('s');  ylabel('\phi');
legend(sprintf('本征值为%.5f的波函数',en));
title('i). \phi(s)的函数曲线');
grid on;
%-------函数定义--------
function y=Numerov(n,h,y0,delta,en)
    y(1)=y0; y(2)=h*delta; 
    con=2*(1-5*h^2*en/6)/(1+h^2*en/6);
    for i=2:n
        y(i+1)=con*y(i)-y(i-1);
    end
end