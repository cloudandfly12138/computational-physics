%pendulum
%2018/10/31 林祥
%单摆(阻力&驱动力)  (thelta)''+q*(thelta)'+g/l*sin(thelta)=b*cos(w0*t)
%理想单摆 q=b=0; w0任意
clc;  clear; format long;
global q b w0 g l
%%----------理想单摆---------------
%q=0.0; b=0.0; w0=2/3;
%%------ 有阻力及驱动力 ---------
%q=0.5;b=0.9;w0=2/3; 
q=0.5;b=1.15;w0=2/3;
l=9.8;   %单摆长度
g=9.8;   %重力加速度
tspan=[1,100];   n=1000;  %时间范围及 步数
y0=[0;1];  % theta, w初值

%%%%%%%%%%%求解%%%%%%%%
[t1,y1]=eluer(tspan,n,y0);
[t2,y2]=RK4(@dfun,tspan,n,y0);  
%%%%%%%%%%%画图%%%%%%%%%%%    
figure(1);   
subplot(2,1,1); set(gca,'Fontsize',16);
plot(t1,y1(:,1),'r.');     hold on;
plot(t2,y2(:,1),'g.');     hold on;
xlabel('t'); ylabel('\theta');
legend('改进欧拉法','RK4法');
title('q=0.5;b=1.15;w0=2/3; \theta随时间变化曲线');
grid on;
subplot(2,1,2);  set(gca,'Fontsize',16);
plot(t1,y1(:,2),'r.');     hold on;
plot(t2,y2(:,2),'g.');     hold on;
xlabel('t'); ylabel('w');
legend('改进欧拉法','RK4法');
title('q=0.5;b=1.15;w0=2/3; 角速度w随时间变化曲线');
grid on;
figure(2);  set(gca,'Fontsize',16);
plot(t1,1-cos(y1(:,1))+(y1(:,2)).^2/2,'r.');  hold on; %机械能是否守恒
plot(t2,1-cos(y2(:,1))+(y2(:,2)).^2/2,'g.');  hold on;%机械能是否守恒
xlabel('t'); ylabel('1-cos\theta+w^2/2');
legend('改进欧拉法','RK4法');
title('机械能是否守恒');
grid on;
%%%%%%%%%%函数定义%%%%%%%%%%%
function [tout,yout]=eluer(tspan,n,y0)
    global q b w0 g l
    t0=tspan(1); tn=tspan(2);  h=(tn-t0)/n;
    t=t0:h:tn;      tout=t';  %列向量
    y0=y0';         yout(1,:)=y0; p(1,:)=y0;  %m(=2)列n+1行
	for i=1:n
	    p(i+1,1)=yout(i,1)+h*yout(i,2);     %显式格式，并作为预测 
        p(i+1,2)=yout(i,2)+h*(b*cos(w0*t(i))-q*yout(i,2)-g/l*sin(yout(i,1)));
        yout(i+1,1)=yout(i,1)+h/2*(yout(i,2)+p(i+1,2));     %用梯形公式校正
        yout(i+1,2)=yout(i,2)+h/2*(b*cos(w0*t(i))+b*cos(w0*t(i+1))-q*yout(i,2)-q*p(i+1,2)-g/l*sin(yout(i,1))-g/l*sin(p(i+1,1)));
    end
end
function [tout,yout]=RK4(fun,tspan,n,y0)
    t0=tspan(1); tn=tspan(2);  h=(tn-t0)/n;
    t=t0;  y=y0(:);   %y为当前步下的函数值(m维列向量)(m为方程组维数)
    tout=t;   yout=y.';
    while(t<tn)
        if t+h>tn, h=tn-t; end
        s1=feval(fun,t,y); s1=s1(:);
        s2=feval(fun,t+h/2,y+h*s1/2);  s2=s2(:);
        s3=feval(fun,t+h/2,y+h*s2/2);  s3=s3(:);
        s4=feval(fun,t+h,y+h*s3);  s4=s4(:);
        t=t+h;
        y=y+h*(s1+2*s2+2*s3+s4)/6;
        tout=[tout;t];
        yout=[yout;y.'];   %m(=2)列n+1行
    end 
end
function dy=dfun(t,y)
    global q b w0 g l
    dy=[y(2) ; b*cos(w0*t)-q*y(2)-g/l*sin(y(1))];
end