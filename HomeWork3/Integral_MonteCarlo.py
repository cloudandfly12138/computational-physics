# Integral_M0nteCarlo
# 2018/12/23 林祥
# 运用平均值法AND投点法计算int(exp(x),0,1)
from random import random
from math import e, exp ,sqrt

x0 = 0;  xn = 1; h = 0.01; N = int((xn-x0)/h)

integral = 0
for i in range(N):
    x = random()*(xn-x0)+x0
    integral += exp(x)
integral /= N
print("平均值法计算int(exp(x),0,1)= %f" % integral)

num = 0
for i in range(N):
    x = random()*(xn-x0)+x0
    y = random()*(exp(xn)-0)
    if y <= exp(x):
        num += 1
integral = num/N * (xn-x0) * (exp(xn)-0)
print("投点法计算int(exp(x),0,1)= %f" % integral)

intergral = 0
for i in range(N):
	x1 = random()*(xn-x0)+x0
	x2 = (-1+ sqrt(1+8*x1) )/2
	y = 2*exp(x2) /(2*x2+1)
	intergral += y
intergral /= N
print("改进平均值法计算int(exp(x),0,1) = %f"% intergral)

print("int(exp(x),0,1)精确解为 %f" % (exp(xn)-exp(x0)))
