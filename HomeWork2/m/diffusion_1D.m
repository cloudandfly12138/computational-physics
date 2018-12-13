%diffusion_1D
%2018/12/4 ����
%u(i,l+1)=u(i,l)+lambda*(u(i+1,l)-2*u(i,l)+u(i-1,l)); ��ʽ
%-lambda*u(i-1,l+1)+(1+2*lambda)*u(i,l+1)-lambda*u(i+1,l+1)=u(i,l); ��ʽ
%u(0,t)=sin(t);  u(1,t)=0; �߽�����
%u(x,0)=0;  ��ʼ����
clc; clear; format long;
t0=0;  tn=1;  tau=0.0001; n=(tn-t0)/tau;
tau_plot=0.01;  n_plot=(tn-t0)/tau_plot;  %��ͼ����ʱ�䲽��
x0=0; xn=1; h=0.02;  N=(xn-x0)/h;
lambda=tau/(h*h); %��������lambda<=1/2;
t_plot=t0:tau_plot:tn;        x_plot=x0:h:xn;
%��ʼ����
u=zeros(N+1,1);  %��ʼ�¶�Ϊ��
u_plot(:,1)=u(:);   iplot=2;

% %��ʽFTCS���
% for l=1:n
%     u(1)=sin(l*tau); u(N+1)=0; %�߽�����
%     u(2:N)=u(2:N)+lambda*(u(3:N+1)-2*u(2:N)+u(1:N-1));
%     %��ͼ��������
%     if(mod(l,n/n_plot)==0)
%         u_plot(:,iplot)=u(:);
%         iplot=iplot+1;
%     end
% end
% figure(1); set(gca,'Fontsize',16);
% mesh(t_plot,x_plot,u_plot);
% xlabel('t');ylabel('x');
% title('һά��ɢ����(��ʽ��)');

%��ʽFTCS��� ���ԽǾ���׷�Ϸ�
for l=1:n
    %ȷ�����ԽǾ���ϵ��
    a=zeros(N+1,1)-lambda;  a(N+1)=0;   %a(1)δʹ��
    b=zeros(N+1,1)+1+2*lambda;  b(1)=1;  b(N+1)=1;
    c=zeros(N+1,1)-lambda;  c(1)=0;   %c(N+1)δʹ��
    u1=tri(a,b,c,u);    u=u1;
    u(1)=sin(l*tau); u(N+1)=0; %�߽���������
    %��ͼ��������
    if(mod(l,n/n_plot)==0)
        u_plot(:,iplot)=u(:);
        iplot=iplot+1;
    end
end
figure(2); set(gca,'Fontsize',16);
mesh(t_plot,x_plot,u_plot);
xlabel('t');ylabel('x');
title('һά��ɢ����(��ʽ��)');

%���ԽǾ�����⺯��
function x=tri(a,b,c,f)
    n=length(f);
    x=zeros(n,1);   y=zeros(n,1);
    d=zeros(n,1);   u=zeros(n,1);
    d(1)=b(1);
    for i=1:n-1
        u(i)=c(i)/d(i);
        d(i+1)=b(i+1)-a(i+1)*u(i);
    end
    %׷
    y(1)=f(1)/d(1);
    for i=2:n
        y(i)=(f(i)-a(i)*y(i-1))/d(i);
    end
    %��
    x(n)=y(n);
    for i=n-1:-1:1
        x(i)=y(i)-u(i)*x(i+1);
    end
end