#diffusion_1D
#2018/12/12 林祥
#u(i,l+1)=u(i,l)+lambda*(u(i+1,l)-2*u(i,l)+u(i-1,l)); 显式
#-lambda*u(i-1,l+1)+(1+2*lambda)*u(i,l+1)-lambda*u(i+1,l+1)=u(i,l); 隐式
#u(0,t)=sin(t);  u(1,t)=0; 边界条件
#u(x,0)=0;  初始条件

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号
#有中文出现的情况，需要u'内容'

def main():
    t0=0.;  tn=1.;  tau=0.0001; n=int((tn-t0)/tau)
    tau_plot=0.01;  n_plot=int((tn-t0)/tau_plot)  #画图所用时间步长
    x0=0.; xn=1.; h=0.02;  N=int((xn-x0)/h)
    lamda=tau/(h*h) #lamda==1/4<=1/2
    t_plot=np.linspace(t0,tn,n_plot+1)
    x_plot=np.linspace(x0,xn,N+1)

    #初始条件
    u=np.zeros(N+1) #初始温度为0
    u_plot=np.zeros((n_plot+1,N+1))
    u_plot[0][:]=u[:];   iplot=1

    # #显式FTCS求解
    # for l in range(1,n+1):
    #     u[0]=np.sin(l*tau); u[-1]=0 #边界条件
    #     u[1:-1]=u[1:-1]+lamda*(u[2:]-2*u[1:-1]+u[0:-2])
    #     #画图所用数据
    #     if(np.mod(l,n/n_plot)==0):
    #         u_plot[iplot][:]=u[:]
    #         iplot+=1
    # fig=plt.figure();   ax=Axes3D(fig)
    # x_plot,t_plot = np.meshgrid(x_plot,t_plot)
    # ax.plot_surface(x_plot,t_plot,u_plot,rstride=1, cstride=1, cmap='rainbow')
    # plt.xlabel('x');    plt.ylabel('t')
    # plt.axis([1,0,0,1])
    # plt.title(u'一维扩散问题(显式法)')
    # plt.show()

    #隐式FTCS求解 三对角矩阵追赶法
    for l in range(1,n+1):
        #确定三对角矩阵系数
        a=np.zeros(N+1)-lamda;  a[-1]=0  #a[0]未使用
        b=np.zeros(N+1)+1+2*lamda;  b[0]=1;  b[-1]=1
        c=np.zeros(N+1)-lamda;  c[0]=0   #c[N]未使用
        u1=tri(a,b,c,u);    u=u1
        u[0]=np.sin(l*tau); u[-1]=0 #边界条件修正
        #画图所用数据
        if(np.mod(l,n/n_plot)==0):
            u_plot[iplot][:]=u[:]
            iplot+=1
    fig=plt.figure();   ax=Axes3D(fig)
    x_plot,t_plot = np.meshgrid(x_plot,t_plot)
    ax.plot_surface(x_plot,t_plot,u_plot,rstride=1, cstride=1, cmap='rainbow')
    plt.xlabel('x');    plt.ylabel('t')
    plt.axis([1,0,0,1])
    plt.title(u'一维扩散问题(隐式法)')
    plt.show()

#三对角矩阵求解函数
def tri(a,b,c,f):
    n=len(f)
    x=np.zeros(n);   y=np.zeros(n)
    d=np.zeros(n);   u=np.zeros(n)
    d[0]=b[0]
    for i in range(n-1):
        u[i]=c[i]/d[i]
        d[i+1]=b[i+1]-a[i+1]*u[i]
    #追
    y[0]=f[0]/d[0]
    for i in range(1,n):
        y[i]=(f[i]-a[i]*y[i-1])/d[i]
    #赶
    x[-1]=y[-1]
    for i in range(n-2,-1,-1):
        x[i]=y[i]-u[i]*x[i+1]
    return x

main()