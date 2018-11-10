%Schrodinger_1D_2
%2018/11/1 ����
%-1/2*phin(s)''+v(s)phin(s)=en*phin(s)
% i). v(s)=0,0<=s<=1; v(s)=inf.s<0 or s>1;
% ii). v(s)=-(1-s^2/2);   

%ii). ����ֵ����Ϊn-1/2;n=0,1,2,3...nΪ�ڵ���
clc;  clear;  format long;
%%t=s; y=phi;
h=0.01;   %ʱ�䲽��
y0=0;  delta=1;    parity=-1;    %�����
%y0=0.3;  delta=0;     parity=1;  %ż���
en=5;  den=0.0002;  %����ֵ�ĳ�ʼ���ֵ�ͳ�ʼ����
ETol=1E-5;
%----------------���--------------
[n,tn,t,yn]=Numerov2(h,y0,delta,en);  %����
while abs(yn(n+1))>ETol  %ָ��˥������֤
    tn=tn+h; n=n+1; t=[t,tn];
    yn(n+1)=2*(1-5*h^2/12*(2*en+2-(t(n))^2))*yn(n)-(1+h^2/12*(2*en+2-(t(n-1))^2))*yn(n-1);
    yn(n+1)=yn(n+1)/(1+h^2/12*(2*en+2-(t(n+1))^2));
    if abs(yn(n+1))>abs(yn(n))   %����˥�����ı�en,���¼���
        en=en+den;     [n,tn,t,yn]=Numerov2(h,y0,delta,en);  
    end
end
%---------------��һ��---------------
sum2=0;  
for i=1:n          %���λ��ֹ�ʽ
    sum2=sum2+h*((yn(i)+yn(i+1))/2)^2;  
end
sum2=2*sum2;   %���ϸ���
yn=yn./sqrt(sum2); 
%----------------��ͼ-----------------
figure(1);  set(gca,'Fontsize',16);
plot(t,yn,'r.',-1.*t,parity.*yn,'r.');   hold on; %��ȫ����
xlabel('s');  ylabel('\phi');
legend(sprintf('����ֵΪ%.4f�Ĳ�����',en));
title('ii). \phi(s)�ĺ�������');
grid on;
%---------��������----------
function [n,tn,t,y]=Numerov2(h,y0,delta,en)%�������������
    y(1)=y0; y(2)=y0+h*delta; 
    tn=sqrt(2*en+2);    n=floor(tn/h);  t=0:h:tn;%����
    for i=2:n
        y(i+1)=2*(1-5*h^2/12*(2*en+2-(t(i))^2))*y(i)-(1+h^2/12*(2*en+2-(t(i-1))^2))*y(i-1);
        y(i+1)=y(i+1)/(1+h^2/12*(2*en+2-(t(i+1))^2));
    end
end