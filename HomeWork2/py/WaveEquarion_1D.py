#WaveEquation_1D
#2018/12/12 林祥
#y(i,l+1)=-y(i,l-1)+(2-2*lambda^2)*y(i,l)+lambda^2*(y(i+1,l)+y(i-1,l));

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号
#有中文出现的情况，需要u'内容'

t0=0;  tn=1;  tau=1E-5; n=int((tn-t0)/tau)
x0=0; xn=1; h=0.01;  N=int((xn-x0)/h);  x=np.linspace(x0,xn,N+1)
v1=300; v2=150
lambda1=v1*tau/h;   lambda2=v2*tau/h

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
for l in range(100):  #时间 求1E-3时range内为100,1E-2时为1000
    for i in range(1,N):
        if(x[i]>0 and x[i]<0.5):
            y[2][i] = (2-2*lambda1**2)*y[1][i] + lambda1**2*(y[1][i+1]+y[1][i-1]) - y[0][i]
        else:
            y[2][i] = (2-2*lambda2**2)*y[1][i] + lambda2**2*(y[1][i+1]+y[1][i-1]) - y[0][i]
    y[2][-1]=y[2][-2] #x=1边界条件
    plt.pause(0.01); plt.clf()  #暂停+清除当前图像
    plt.plot(x,y[2]) 
    plt.xlabel('x'); plt.ylabel('y'); plt.grid()
    plt.axis([0,1,-0.5,0.5])
    t=(l+1)*tau;  plt.title(u't=%.5f时刻波形图'%t)
    y[0]=y[1];  y[1]=y[2]
plt.ioff() #plot关闭交互模式
plt.show() #显示最后一幅图