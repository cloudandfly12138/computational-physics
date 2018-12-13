%WaveEquation_1D
%2018/12/4 林祥
%y(i,l+1)=-y(i,l-1)+(2-2*lambda^2)*y(i,l)+lambda^2*(y(i+1,l)+y(i-1,l));
clc; clear; format long;
t0=0;  tn=1;  tau=1E-5; n=(tn-t0)/tau;
x0=0; xn=1; h=0.01;  N=(xn-x0)/h;  x=x0:h:xn;
v1=300; v2=150;
lambda1=v1*tau/h;   lambda2=v2*tau/h;
%初始条件 和 边界条件
y=zeros(N+1,3);
for i=1:N+1
    y(i,1)=exp(-1000*(x(i)-0.3)^2);
    y(i,2)=y(i,1);
end
y(1,:)=0;  

%迭代求解
for l=1:100 %时间迭代
    for i=2:N
        if x(i)>0 && x(i)<0.5
            y(i,3)=-y(i,1)+(2-2*lambda1^2)*y(i,2)+lambda1^2*(y(i+1,2)+y(i-1,2));
        else
            y(i,3)=-y(i,1)+(2-2*lambda2^2)*y(i,2)+lambda2^2*(y(i+1,2)+y(i-1,2));
        end
    end
    y(N+1,3)=y(N,3);%边界条件
    figure(2); set(gca,'Fontsize',16);
    plot(x,y(:,3));
    xlabel('x');ylabel('y'); grid on;
    axis([0,1,-0.5,0.5]);
    title(sprintf('t=%.5f时刻波形图',l*tau));
    y(:,1)=y(:,2);  y(:,2)=y(:,3); 
end