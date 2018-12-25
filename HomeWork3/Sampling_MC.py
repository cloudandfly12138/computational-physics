# Sampling_MC
# 2018/12/25 林祥
# 对f(x)=sqrt(2/pi)*exp(-x**2/2) (0<=x<inf) 进行抽样
from random import random
from math import e, pi, exp, sqrt, log

N = 100
# 舍选法 反函数法结合
y = []
for i in range(N):
    x1 = random()
    x2 = -log(x1,e)
    x3 = random()
    if x3 <= sqrt(2*e/pi)*exp(-(x2-1)**2/2):
        y.append(x2)
print("抽样结果:",end=' ')
print(",".join(str(i) for i in y))
