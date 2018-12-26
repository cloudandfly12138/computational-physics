#WaveEquation_1D
#2018/12/12 林祥
#y(i,l+1)=-y(i,l-1)+(2-2*lambda^2)*y(i,l)+lambda^2*(y(i+1,l)+y(i-1,l));

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号
#有中文出现的情况，需要u'内容'

t0=0;  tn=1;  tau=1.25E-5; n=int((tn-t0)/tau)
x0=0; xi=0.5; xn=1; h1=0.005; h2=0.0025  #x在(0,0.5)和(0.5,1)采用不同步长
N1=int((xi-x0)/h1); N2=int((xn-xi)/h2); N=N1+N2
x1=np.linspace(x0,xi,N1+1); x2=np.linspace(xi+h2,xn,N2)
x=np.r_[x1,x2]
v1=300; v2=150
lambda1=v1*tau/h1;   lambda2=v2*tau/h2 #lambda1==lambda2==0.75 在（sqrt(2)/2, 1）之间

#初始条件 和 边界条件
y=np.zeros((3,N+1))  #初始化
for i in range(N+1):
    y[0][i]=np.exp(-1000*(x[i]-0.3)**2)
    y[1][i]=y[0][i]
y[:,0]=0   #x=0边界条件 #de了我一个晚上将加半个下午的BUG!!!!!
#也可写为y[0][0],y[1][0],y[2][0]=0,0,0   #不能写为y[:][0]=0,会把第0行(非列)清零

fig=plt.figure()
plt.ion() #plot打开交互模式
#求解
for l in range(800):  #时间 求t=1E-3时range内为80,1E-2时为800
    for i in range(1,N):
        if(x[i]>0 and x[i]<xi):  #lambda1==lambda2,实际不用分支结构就行.
            y[2][i] = (2-2*lambda1**2)*y[1][i] + lambda1**2*(y[1][i+1]+y[1][i-1]) - y[0][i]
        else:
            y[2][i] = (2-2*lambda2**2)*y[1][i] + lambda2**2*(y[1][i+1]+y[1][i-1]) - y[0][i]
    y[2][-1]=y[2][-2] #x=1边界条件
    plt.pause(0.01); plt.clf()  #暂停+清除当前图像
    plt.plot(x,y[2]) 
    plt.xlabel('x'); plt.ylabel('y'); plt.grid()
    plt.axis([0,1,-1,1])
    t=(l+1)*tau;  plt.title(u't=%.5f时刻波形图'%t)
    y[0]=y[1];  y[1]=y[2]
plt.ioff() #plot关闭交互模式
plt.show() #显示最后一幅图
