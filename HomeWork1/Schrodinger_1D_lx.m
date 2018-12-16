%Schrodinger_1D
%2018/11/1 ����
%-1/2*phin(s)''+v(s)phin(s)=en*phin(s)
% i). v(s)=0,0<=s<=1; v(s)=inf.s<0 or s>1;
% ii). v(s)=-(1-s^2/2);
clc;  clear;  format long;
%%t=s; y=phi;
h=0.01;  n=1/h;   t=0:h:1;  %ʱ�䷶Χ�� ����
y0=0; yn=0;  %��ֵ
delta=1;  %��һ����������!?, delta��Ϊ����ֵ
en=20;  den=0.01;  %����ֵ�ĳ�ʼ���ֵ�ͳ�ʼ����
ETol=1E-5;
%----------------���--------------
yn1=Numerov(n,h,y0,delta,en);  %tspan=[0,1]
while   abs(den)>ETol
    en=en+den;
    yn2=Numerov(n,h,y0,delta,en);
    if   (yn2(n+1))*(yn1(n+1))>0  %��en������ʵ��yn=0�ı�ֵ����
    else
        en=en-den;
        den=den/2;
    end
end
%---------------��һ��---------------
sum1=0;  
for i=1:n          %���λ��ֹ�ʽ
    sum1=sum1+h*((yn2(i)+yn2(i+1))/2)^2;  
end
yn2=yn2./sqrt(sum1);    
%----------------��ͼ-----------------
figure(1);  set(gca,'Fontsize',16);
plot(t,yn2,'r.');   hold on;
xlabel('s');  ylabel('\phi');
legend(sprintf('����ֵΪ%.5f�Ĳ�����',en));
title('i). \phi(s)�ĺ�������');
grid on;
%-------��������--------
function y=Numerov(n,h,y0,delta,en)
    y(1)=y0; y(2)=h*delta; 
    con=2*(1-5*h^2*en/6)/(1+h^2*en/6);
    for i=2:n
        y(i+1)=con*y(i)-y(i-1);
    end
end