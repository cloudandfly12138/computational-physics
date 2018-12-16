%cannonshell
%2018/10/30 林祥
clc;  clear; format long;
issue=3;            %是否考虑阻力，no[1]，yes[2],both[0],寻找击中(15k,3k)的发射角[3]
t0=0;   tn=100;   dt=1;  n=(tn-t0)/dt;
t=t0:dt:tn;
v0=500; %初始速度
g=9.8;   %重力加速度    
b=2E-5;             %B2/m
if(issue==3)
    thelta=atan(1/5);     %issue3的初始发射角
else
    thelta=60;  %初始发射角(以度为单位)
    thelta=thelta*pi/180;  
end
%y1=x, y2=z, y3=vx, y4=vz;
y10=0; y20=0; y30=v0*cos(thelta); y40=v0*sin(thelta);  %初值

%%%%%%%%%%%求解%%%%%%%%
[y11,y21,y31,y41]=eluer(n,dt,y10,y20,y30,y40,v0,0); %无阻力
[y12,y22,y32,y42]=eluer(n,dt,y10,y20,y30,y40,v0,b);  %有阻力

if(issue==3)
    dthelta=0.02; %寻找合适发射角时thelta的变化量0.02*180/pi=1.146度
    ETol=0.01*pi/180;  %寻找到的发射角最小误差不超过0.01度
    k=0;        %迭代次数
    while(abs(dthelta)>ETol && thelta<pi/2 && thelta>0 )
        k=k+1;
        thelta=thelta+dthelta;
        y10=0; y20=0; y30=v0*cos(thelta); y40=v0*sin(thelta);  %初值
        [y12,y22,y32,y42]=eluer(n,dt,y10,y20,y30,y40,v0,b);
        figure(1);
        plot(y12,y22,'g.');     
        legend('有阻力(改进欧拉法)');
        for i=1:n
            if((y12(i)-15000)*(y12(i+1)-15000)<0)  %找到x=15000的点
                break;
            end
        end
        if((y12(i)+y12(i+1))/2<15000)   i=i+1;    %找到x=15000的点
        end
        if(y22(i)>3000)        
            thelta=thelta-dthelta;
            dthelta=dthelta/2;
        end
    end
   sprintf('有阻力下要击中(15k,3k) ,发射角为%0.2f度，误差不超过0.01度,迭代次数为%d',thelta*180/pi,k)
end

%%%%%%%%%%%画图%%%%%%%%%%%    
figure(1);
set(gca,'Fontsize',16);
if(issue==0)
    plot(y11,y21,'r.');     hold on;
    plot(v0*cos(thelta).*t,v0*t*sin(thelta)-g/2.*t.^2,'b-');        hold on;
    plot(y12,y22,'g.');     hold on;
    legend('无阻力(改进欧拉法)','无阻力(解析解)','有阻力(改进欧拉法)');
elseif(issue==1)
    plot(y11,y21,'r.');     hold on;
    plot(v0*cos(thelta).*t,v0*t*sin(thelta)-g/2.*t.^2,'b-');        hold on;
    legend('无阻力(改进欧拉法)','无阻力(解析解)');
elseif(issue==2)
    plot(y12,y22,'g.');     hold on;
     legend('有阻力(改进欧拉法)');
end
xlabel('x'); ylabel('z');
title('cannonshell炮弹轨迹');
grid on;

%%%%%%%%%%函数定义%%%%%%%%%%%
function [y1,y2,y3,y4]=eluer(n,dt,y10,y20,y30,y40,v0,b)    %有阻力下改进欧拉方法
        p1(1)=y10;  p2(1)=y20;  p3(1)=y30;  p4(1)=y40;    v=v0;
        y1(1)=y10;   y2(1)=y20;  y3(1)=y30;  y4(1)=y40;
        g=9.8;   %重力加速度
        for i=1:n                      %有阻力
            p1(i+1)=y1(i)+dt*y3(i);     %显式格式，并作为预测
            p2(i+1)=y2(i)+dt*y4(i);
            p3(i+1)=y3(i)-dt*b*v*y3(i);
            p4(i+1)=y4(i)-dt*(g+b*v*y4(i));
            y1(i+1)=y1(i)+dt/2*(y3(i)+p3(i+1));     %用梯形公式校正
            y2(i+1)=y2(i)+dt/2*(y4(i)+p4(i+1));
            y3(i+1)=y3(i)-dt/2*(b*v*y3(i)+b*v*p3(i+1));
            y4(i+1)=y4(i)-dt/2*(2*g+b*v*y4(i)+b*v*p4(i+1));
            v=sqrt((y3(i+1))^2+(y4(i+1))^2);  %新的速度
        end
end