%pendulum
%2018/10/31 ����
%����(����&������)  (thelta)''+q*(thelta)'+g/l*sin(thelta)=b*cos(w0*t)
%���뵥�� q=b=0; w0����
clc;  clear; format long;
global q b w0 g l
%%----------���뵥��---------------
%q=0.0; b=0.0; w0=2/3;
%%------ �������������� ---------
%q=0.5;b=0.9;w0=2/3; 
q=0.5;b=1.15;w0=2/3;
l=9.8;   %���ڳ���
g=9.8;   %�������ٶ�
tspan=[1,100];   n=1000;  %ʱ�䷶Χ�� ����
y0=[0;1];  % theta, w��ֵ

%%%%%%%%%%%���%%%%%%%%
[t1,y1]=eluer(tspan,n,y0);
[t2,y2]=RK4(@dfun,tspan,n,y0);  
%%%%%%%%%%%��ͼ%%%%%%%%%%%    
figure(1);   
subplot(2,1,1); set(gca,'Fontsize',16);
plot(t1,y1(:,1),'r.');     hold on;
plot(t2,y2(:,1),'g.');     hold on;
xlabel('t'); ylabel('\theta');
legend('�Ľ�ŷ����','RK4��');
title('q=0.5;b=1.15;w0=2/3; \theta��ʱ��仯����');
grid on;
subplot(2,1,2);  set(gca,'Fontsize',16);
plot(t1,y1(:,2),'r.');     hold on;
plot(t2,y2(:,2),'g.');     hold on;
xlabel('t'); ylabel('w');
legend('�Ľ�ŷ����','RK4��');
title('q=0.5;b=1.15;w0=2/3; ���ٶ�w��ʱ��仯����');
grid on;
figure(2);  set(gca,'Fontsize',16);
plot(t1,1-cos(y1(:,1))+(y1(:,2)).^2/2,'r.');  hold on; %��е���Ƿ��غ�
plot(t2,1-cos(y2(:,1))+(y2(:,2)).^2/2,'g.');  hold on;%��е���Ƿ��غ�
xlabel('t'); ylabel('1-cos\theta+w^2/2');
legend('�Ľ�ŷ����','RK4��');
title('��е���Ƿ��غ�');
grid on;
%%%%%%%%%%��������%%%%%%%%%%%
function [tout,yout]=eluer(tspan,n,y0)
    global q b w0 g l
    t0=tspan(1); tn=tspan(2);  h=(tn-t0)/n;
    t=t0:h:tn;      tout=t';  %������
    y0=y0';         yout(1,:)=y0; p(1,:)=y0;  %m(=2)��n+1��
	for i=1:n
	    p(i+1,1)=yout(i,1)+h*yout(i,2);     %��ʽ��ʽ������ΪԤ�� 
        p(i+1,2)=yout(i,2)+h*(b*cos(w0*t(i))-q*yout(i,2)-g/l*sin(yout(i,1)));
        yout(i+1,1)=yout(i,1)+h/2*(yout(i,2)+p(i+1,2));     %�����ι�ʽУ��
        yout(i+1,2)=yout(i,2)+h/2*(b*cos(w0*t(i))+b*cos(w0*t(i+1))-q*yout(i,2)-q*p(i+1,2)-g/l*sin(yout(i,1))-g/l*sin(p(i+1,1)));
    end
end
function [tout,yout]=RK4(fun,tspan,n,y0)
    t0=tspan(1); tn=tspan(2);  h=(tn-t0)/n;
    t=t0;  y=y0(:);   %yΪ��ǰ���µĺ���ֵ(mά������)(mΪ������ά��)
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
        yout=[yout;y.'];   %m(=2)��n+1��
    end 
end
function dy=dfun(t,y)
    global q b w0 g l
    dy=[y(2) ; b*cos(w0*t)-q*y(2)-g/l*sin(y(1))];
end