# laplace_2D
# 2018/12/12 林祥
# u(i,j)=(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1))/4;
# issue==1 u(0,y)=0; u(1,y)=0; u(x,0)=0; u(x,1)=100; 
# issue==2 u(0,y)=0; u(1,y)=0; u(x,0)=0; u(x,1)=0;  u(0.3:0.7,0.4)=-100; u(0.3:0.7,0.6)=100;

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
plt.rcParams['axes.unicode_minus']=False #用来正常显示负号
#有中文出现的情况，需要u'内容'

def main():
    x0=0; xn=1; y0=0; yn=1; h=0.01
    Mx=int((xn-x0)/h); My=int((yn-y0)/h)
    tol=1E-5; MaxIter=10000  #容忍误差 最大迭代次数
    lamda=1.5  #超松弛迭代因子
    u=np.zeros((My+1,Mx+1))  #初始化
    u0=np.zeros((My+1,Mx+1))  #初始化u0

    # #边界条件1
    # issue=1
    # for j in range(My+1):
    #     u[j][0],u[j][-1]=0,0
    # for i in range(Mx+1):
    #     u[0][i],u[-1][i]=0,100
    # #边界平均值作迭代初值1
    # sum_of_bv = np.sum(u[1:-1][0])+np.sum(u[1:-1][-1])+np.sum(u[0][1:-1])+np.sum(u[-1][1:-1]) 
    # average=sum_of_bv/(2*(Mx+My-2))
    # for j in range(1,My):
    #     for i in range(1,Mx):
    #         u[j][i]=average

    #边界条件2 迭代初值2
    issue=2
    u=np.zeros((My+1,Mx+1)) #迭代初值2
    for j in range(My+1):
        u[j][0],u[j][-1]=0,0
    for i in range(Mx+1):
        u[0][i],u[-1][i]=0,0
    for i in range( int((0.3-x0)/h) , int((0.7-x0)/h+1) ):
        u[int((0.4-y0)/h)][i] , u[int((0.6-y0)/h)][i] = -100,100

    #迭代
    for itr in range(MaxIter):
        for i in range(1,My):
            for j in range(1,Mx):
                if(issue==1 or issue==2 and not(j>=(0.3-x0)/h and j<=(0.7-x0)/h and (i==(0.4-y0)/h or i==(0.6-y0)/h) ) ):
                    u[i][j]=lamda*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1])/4+(1-lamda)*u[i][j] #超松弛迭代
                    #u[i][j]=(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1])/4
        if (itr>0 and np.max(np.abs(u-u0))<tol ):  
            break
        u0=u.copy() #不能直接赋值!!!!!!!
    print('itr=%d,tol=%.8f'%(itr,tol))
    #画图
    fig=plt.figure();   ax=Axes3D(fig)
    x=np.linspace(x0,xn,Mx+1); y=np.linspace(y0,yn,My+1)
    X,Y = np.meshgrid(x,y)
    ax.plot_surface(X,Y,u,rstride=1, cstride=1, cmap='rainbow')
    plt.xlabel('x');  plt.ylabel('y')
    plt.axis([1,0,0,1])
    plt.title(u'二维静电场问题(laplace)(问题%d)'%issue)
    plt.show()

main()