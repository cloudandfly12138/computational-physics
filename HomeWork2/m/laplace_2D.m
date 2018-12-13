%laplace_2D
%2018/12/4 ����
%u(i,j)=(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1))/4;
%issue==1 u(0,y)=0; u(1,y)=0; u(x,0)=0; u(x,1)=100; 
%issue==2 u(0,y)=0; u(1,y)=0; u(x,0)=0; u(x,1)=0;  u(0.3:0.7,0.4)=-100; u(0.3:0.7,0.6)=100;
clc; clear; format long;
x0=0; xn=1; y0=0; yn=1; h=0.01;
Mx=(xn-x0)/h; My=(yn-y0)/h;
tol=1E-5; MaxIter=10000; %������� ����������
lambda=1.5;  %���ɳڵ�������

% %�߽�����1
% issue=1;
% for j=1:My+1
%     u(j,[1 Mx+1])=[0 0];
% end
% for i=1:Mx+1
%     u([1 My+1],i)=[0 ;100];
% end
% %�߽�ƽ��ֵ��������ֵ1
% sum_of_bv=sum(sum([u(2:My,[1 Mx+1])   u([1 My+1],2:Mx)']));
% u(2:My,2:Mx)=sum_of_bv/(2*(Mx+My-2));

%�߽�����2 ������ֵ2
issue=2;
u=zeros(My+1,Mx+1);  %������ֵ2
for j=1:My+1
    u(j,[1 Mx+1])=[0 0];
end
for i=1:Mx+1
    u([1 My+1],i)=[0; 0];
end
for i=(0.3-x0)/h+1 : (0.7-x0)/h+1
    u([(0.4-y0)/h+1   (0.6-y0)/h+1],i)=[-100 ;100];
end

%����
for itr=1:MaxIter
    for i=2:My
        for j=2:Mx
            if issue==1 || ( issue==2 && ~(j>=(0.3-x0)/h+1 && j<=(0.7-x0)/h+1 && (i==(0.4-y0)/h+1 ||i==(0.6-y0)/h+1) ) )
                u(i,j)=lambda*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1))/4+(1-lambda)*u(i,j);%���ɳڵ���
                %u(i,j)=(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1))/4;
            end
        end
    end
    if itr>1 && max(max(abs(u-u0)))<tol
        break
    end
    u0=u;
end
%��ͼ
figure(issue); set(gca,'Fontsize',16);
[X,Y]=meshgrid(x0:h:xn,y0:h:yn);
mesh(X,Y,u);
xlabel('x');ylabel('y');zlabel('\phi');
title(sprintf('��ά���糡����(laplace)(����%d)',issue));

