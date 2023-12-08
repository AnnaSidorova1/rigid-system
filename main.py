import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import math

#система Ван-дер-Поля
# x' = y
# y' = mu*(1 - x*x)*y - x

def f_x(x, y):
    return y

def f_y(x, y):
    return mu*(1-x*x)*y - x

n = 30000
T = 30
x_RK = [None] * (n) 
y_RK = [None] * (n) 

x_G = [None] * (n) 
y_G = [None] * (n) 

h = T / n 
mu = 3
eps = 0.0001

print('h =', h)
x_RK[0] = 2
y_RK[0] = 0

x_G[0] = 2
y_G[0] = 0

#метод Рунге-Кутта    
for i in range(0, n-1):    
    K1 = h*(y_RK[i])
    K2 = h*(y_RK[i] + K1*0.5)
    K3 = h*(y_RK[i] + K2*0.5)
    K4 = h*(y_RK[i] + K3)
    x_RK[i+1] = x_RK[i] + K1 / 6 + K2 / 3 + K3 / 3 + K4 / 6
    
    K1 = h*(mu*(1-x_RK[i]*x_RK[i])*y_RK[i] - x_RK[i])
    K2 = h*(mu*(1-(x_RK[i]+h*0.5)*(x_RK[i]+h*0.5))*(y_RK[i] + K1*0.5) - (x_RK[i]+h*0.5))
    K3 = h*(mu*(1-(x_RK[i]+h*0.5)*(x_RK[i]+h*0.5))*(y_RK[i] + K2*0.5) - (x_RK[i]+h*0.5))
    K4 = h*(mu*(1-(x_RK[i]+h)*(x_RK[i]+h))*(y_RK[i] + K3) - (x_RK[i]+h))
    y_RK[i+1] = y_RK[i] + K1 / 6 + K2 / 3 + K3 / 3 + K4 / 6


#вычислим начальное приближение для метода Гира
x_G[1] = x_RK[1]
y_G[1] = y_RK[1]

x_G[2] = x_RK[2]
y_G[2] = y_RK[2]

x_G[3] = x_RK[3]
y_G[3] = y_RK[3]

#метод Гира
#ищем x_G[i] через приближенное x_G[i], найденное по схеме прогноза
t = 3*h
i = 4
W = [[0] * (2) for _ in range(2)]
R = [[0] * (2) for _ in range(2)]
F = [0] * 2

while(i < n):
    
    #поиск начального приближения
    x_G[i] = (-10/3)*x_G[i-1] + 6*x_G[i-2] - 2*x_G[i-3] + (1/3)*x_G[i-4] + 4*h*f_x(x_G[i-1], y_G[i-1])
    y_G[i] = (-10/3)*y_G[i-1] + 6*y_G[i-2] - 2*y_G[i-3] + (1/3)*y_G[i-4] + 4*h*f_y(x_G[i-1], y_G[i-1])
    
    while True:
        #поиск обратной матрицы и значения функций в точках начального приближения
        det_W = 625 - 300*h*mu - 25*mu*x_G[i]*x_G[i] + 12*12*h*h - 12*h*mu*y_G[i]
        W[0][0] = 25-12*h*mu-mu*x_G[i]*x_G[i]
        W[1][0] = mu*y_G[i] - 12*h
        W[0][1] = 12*h
        W[1][1] = 25
        det_W_rev = 1 / det_W
    
        F[0] = 25*x_G[i] - 12*h*y_G[i] - 48*x_G[i-1] + 36*x_G[i-2] - 16*x_G[i-3] + 3*x_G[i-4] 
        F[1] = 25*y_G[i] - 12*h*mu*y_G[i] + 12*h*mu*x_G[i]*x_G[i]*y_G[i] + 12*h*x_G[i] - 48*y_G[i-1] + 36*y_G[i-2] - 16*y_G[i-3] + 3*y_G[i-4]    
    
        #вычисление произведения матрицы на вектор значения функций
        R[0] = (W[0][0]*F[0] + W[0][1]*F[1]) * det_W_rev
        R[1] = (W[1][0]*F[0] + W[1][1]*F[1]) * det_W_rev
    
        #вычисление итогового результата x ~(k+1)
        x_prev = x_G[i]
        y_prev = y_G[i]
        x_G[i] = x_G[i] - R[0]
        y_G[i] = y_G[i] - R[1]

        if (abs(abs(x_G[i]) - abs(x_prev)) < eps and abs(abs(y_G[i]) - abs(y_prev)) < eps):
            break
        
    i += 1      


 
o_x2 = np.arange(0, T, h)
o_y2 = np.array(x_G)
o_y_22 = np.array(x_RK)
plt.plot(o_x2, o_y2, color ="red")
plt.plot(o_x2, o_y_22, ':', color ="black")
plt.title('y1')

plt.show() 

o_x3 = np.arange(0, T, h)
o_y3 = np.array(y_G)
o_y_33 = np.array(y_RK)
plt.plot(o_x3, o_y3, color ="red")
plt.plot(o_x2, o_y_33, ':', color ="black")
plt.title('y2')

plt.show() 

o_x4 = np.array(x_G)
o_y4 = np.array(y_G)
o_x_44 = np.array(x_RK)
o_y_44 = np.array(y_RK)
plt.plot(o_x4, o_y4, color ="blue")
plt.plot(o_x_44, o_y_44, ':', color ="red")
plt.title('y1 - y2')

plt.show() 
    
    
