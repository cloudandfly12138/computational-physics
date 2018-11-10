%cannonshell
%2018/10/30 ����
clc;  clear; format long;
issue=3;            %�Ƿ���������no[1]��yes[2],both[0],Ѱ�һ���(15k,3k)�ķ����[3]
t0=0;   tn=100;   dt=1;  n=(tn-t0)/dt;
t=t0:dt:tn;
v0=500; %��ʼ�ٶ�
g=9.8;   %�������ٶ�    
b=2E-5;             %B2/m
if(issue==3)
    thelta=atan(1/5);     %issue3�ĳ�ʼ�����
else
    thelta=60;  %��ʼ�����(�Զ�Ϊ��λ)
    thelta=thelta*pi/180;  
end
%y1=x, y2=z, y3=vx, y4=vz;
y10=0; y20=0; y30=v0*cos(thelta); y40=v0*sin(thelta);  %��ֵ

%%%%%%%%%%%���%%%%%%%%
[y11,y21,y31,y41]=eluer(n,dt,y10,y20,y30,y40,v0,0); %������
[y12,y22,y32,y42]=eluer(n,dt,y10,y20,y30,y40,v0,b);  %������

if(issue==3)
    dthelta=0.02; %Ѱ�Һ��ʷ����ʱthelta�ı仯��0.02*180/pi=1.146��
    ETol=0.01*pi/180;  %Ѱ�ҵ��ķ������С������0.01��
    k=0;        %��������
    while(abs(dthelta)>ETol && thelta<pi/2 && thelta>0 )
        k=k+1;
        thelta=thelta+dthelta;
        y10=0; y20=0; y30=v0*cos(thelta); y40=v0*sin(thelta);  %��ֵ
        [y12,y22,y32,y42]=eluer(n,dt,y10,y20,y30,y40,v0,b);
        figure(1);
        plot(y12,y22,'g.');     
        legend('������(�Ľ�ŷ����)');
        for i=1:n
            if((y12(i)-15000)*(y12(i+1)-15000)<0)  %�ҵ�x=15000�ĵ�
                break;
            end
        end
        if((y12(i)+y12(i+1))/2<15000)   i=i+1;    %�ҵ�x=15000�ĵ�
        end
        if(y22(i)>3000)        
            thelta=thelta-dthelta;
            dthelta=dthelta/2;
        end
    end
   sprintf('��������Ҫ����(15k,3k) ,�����Ϊ%0.2f�ȣ�������0.01��,��������Ϊ%d',thelta*180/pi,k)
end

%%%%%%%%%%%��ͼ%%%%%%%%%%%    
figure(1);
set(gca,'Fontsize',16);
if(issue==0)
    plot(y11,y21,'r.');     hold on;
    plot(v0*cos(thelta).*t,v0*t*sin(thelta)-g/2.*t.^2,'b-');        hold on;
    plot(y12,y22,'g.');     hold on;
    legend('������(�Ľ�ŷ����)','������(������)','������(�Ľ�ŷ����)');
elseif(issue==1)
    plot(y11,y21,'r.');     hold on;
    plot(v0*cos(thelta).*t,v0*t*sin(thelta)-g/2.*t.^2,'b-');        hold on;
    legend('������(�Ľ�ŷ����)','������(������)');
elseif(issue==2)
    plot(y12,y22,'g.');     hold on;
     legend('������(�Ľ�ŷ����)');
end
xlabel('x'); ylabel('z');
title('cannonshell�ڵ��켣');
grid on;

%%%%%%%%%%��������%%%%%%%%%%%
function [y1,y2,y3,y4]=eluer(n,dt,y10,y20,y30,y40,v0,b)    %�������¸Ľ�ŷ������
        p1(1)=y10;  p2(1)=y20;  p3(1)=y30;  p4(1)=y40;    v=v0;
        y1(1)=y10;   y2(1)=y20;  y3(1)=y30;  y4(1)=y40;
        g=9.8;   %�������ٶ�
        for i=1:n                      %������
            p1(i+1)=y1(i)+dt*y3(i);     %��ʽ��ʽ������ΪԤ��
            p2(i+1)=y2(i)+dt*y4(i);
            p3(i+1)=y3(i)-dt*b*v*y3(i);
            p4(i+1)=y4(i)-dt*(g+b*v*y4(i));
            y1(i+1)=y1(i)+dt/2*(y3(i)+p3(i+1));     %�����ι�ʽУ��
            y2(i+1)=y2(i)+dt/2*(y4(i)+p4(i+1));
            y3(i+1)=y3(i)-dt/2*(b*v*y3(i)+b*v*p3(i+1));
            y4(i+1)=y4(i)-dt/2*(2*g+b*v*y4(i)+b*v*p4(i+1));
            v=sqrt((y3(i+1))^2+(y4(i+1))^2);  %�µ��ٶ�
        end
end