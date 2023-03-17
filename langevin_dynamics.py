import math
import sys
from math import exp
import numpy as np

def potential(x,y,pes_type):
    
    # For PES-1
    if pes_type == 1:
        return 0.02*(x**4+y**4) - 4*exp(-((x+2)**2 + (y+2)**2)) - 4*exp(-((x-2)**2 + (y-2)**2)) + 0.3*(x-y)**2 + 0.0026
    
    # For PES-2 
    if pes_type == 2:
        return 0.03*(x**4+y**4) - 4*exp(-((x+2)**2 + (y+2)**2)) - 4*exp(-((x-2)**2 + (y-2)**2)) + 0.4*(x-y)**2 + 4*exp(-((x)**2+(y)**2)) - 2.1245
    
    # For PES-3
    if pes_type == 3:
        return 0.02*(x**4+y**4) - 3.73*exp(-((x+2)**2/8 + (y+2)**2/8)) - 3.73*exp(-((x-2)**2/8 + (y-2)**2/8)) + 3*exp(-((x)**2/2 + (y)**2/15)) + 2*exp(-((x)**2/2 + (y)**2/2)) - 0.5085
    
    # For PES-4
    if pes_type == 4:

        # Known constants
        A = [-200, -100, -170, 15]
        a = [-1, -1, -6.5, 0.7]
        b = [0, 0, 11, 0.6]
        c = [-10, -10, -6.5, 0.7]
        x0 = [1, 0, -0.5, -1]
        y0 = [0, 0.5, 1.5, 1]

        # Muller-Brown Potential Equation - Kob 2017 arXiv (https://arxiv.org/pdf/1701.01241.pdf)
      
        return A[0]*exp(a[0]*((x-x0[0])**2) + b[0]*(x-x0[0])*(y-y0[0]) + c[0]*((y-y0[0])**2)) + A[1]*exp(a[1]*((x-x0[1])**2) + b[1]*(x-x0[1])*(y-y0[1]) + c[1]*((y-y0[1])**2)) + A[2]*exp(a[2]*((x-x0[2])**2) + b[2]*(x-x0[2])*(y-y0[2]) + c[2]*((y-y0[2])**2)) + A[3]*exp(a[3]*((x-x0[3])**2) + b[3]*(x-x0[3])*(y-y0[3]) + c[3]*((y-y0[3])**2))
    else:
        print('Invalid choice of pes_type. Choose a number between 1 and 4.')

def kinetic_energy(phasepoint):
    px = phasepoint[2]
    py = phasepoint[3]
    return 0.5*(px*px+py*py) 
 
def force(x,y,px,py,dt,beta,gamma,pes_type):
    std_dev = math.sqrt(2.0*gamma/(beta*dt))
    
    #For PES-1
    if pes_type == 1:
        dV_dx = 0.08*x**3 + 8*(x-2)*exp(-(x-2)**2-(y-2)**2) + 8*(x+2)*exp(-(x+2)**2-(y+2)**2) + 0.6*(x-y)
        dV_dy = -0.6*(x-y) + 8*(y-2)*exp(-(x-2)**2-(y-2)**2) + 8*(y+2)*exp(-(x+2)**2-(y+2)**2) + 0.08*y**3
    
    # For PES-3
    if pes_type == 3:
        dV_dx = 0.08*x**3 - 2*x*exp(-x**2/2-y**2/2) - 3*x*exp(-x**2/2-y**2/15) - 0.9325*(-x-2)*exp(-(1/8)*(x+2)**2-(1/8)*(y+2)**2) - 0.9325*(2-x)*exp(-(1/8)*(x-2)**2 - (1/8)*(y-2)**2)
        dV_dy = -2*y*exp(-x**2/2 - y**2/2) - (2/5)*y*exp(-x**2/2 - y**2/15) - 0.9325*(-y-2)*exp(-(1/8)*(x+2)**2 - (1/8)*(y+2)**2) - 0.9325*(2-y)*exp(-(1/8)*(x-2)**2 - (1/8)*(y-2)**2) + 0.08*y**3
    
    # For PES-2
    if pes_type == 2:
        dV_dx = 0.12*x**3 - 8*x*exp(-x**2-y**2) + 8*(x-2)*exp(-(x-2)**2-(y-2)**2) + 8*(x+2)*exp(-(x+2)**2-(y+2)**2) + 0.8*x - 0.8*y
        dV_dy = 0.12*y**3 - 8*y*exp(-x**2-y**2) + 8*(y-2)*exp(-(x-2)**2-(y-2)**2) + 8*(y+2)*exp(-(x+2)**2-(y+2)**2) - 0.8*x + 0.8*y 
    
    # For PES-4 
    if pes_type == 4:
        
        # Known constants
        A = [-200, -100, -170, 15]
        a = [-1, -1, -6.5, 0.7]
        b = [0, 0, 11, 0.6]
        c = [-10, -10, -6.5, 0.7]
        x0 = [1, 0, -0.5, -1]
        y0 = [0, 0.5, 1.5, 1]
        
        dV_dx = A[0]*(a[0]*2*(x-x0[0]) + b[0]*(y-y0[0]))*exp(a[0]*((x-x0[0])**2) + b[0]*(x-x0[0])*(y-y0[0]) + c[0]*((y-y0[0])**2)) + A[1]*(a[1]*2*(x-x0[1]) + b[1]*(y-y0[1]))*exp(a[1]*((x-x0[1])**2) + b[1]*(x-x0[1])*(y-y0[1]) + c[1]*((y-y0[1])**2)) + A[2]*(a[2]*2*(x-x0[2]) + b[2]*(y-y0[2]))*exp(a[2]*((x-x0[2])**2) + b[2]*(x-x0[2])*(y-y0[2]) + c[2]*((y-y0[2])**2)) + A[3]*(a[3]*2*(x-x0[3]) + b[3]*(y-y0[3]))*exp(a[3]*((x-x0[3])**2) + b[3]*(x-x0[3])*(y-y0[3]) + c[3]*((y-y0[3])**2)) + A[3]*(a[3]*2*(x-x0[3]) + b[3]*(y-y0[3]))*exp(a[3]*((x-x0[3])**2) + b[3]*(x-x0[3])*(y-y0[3]) + c[3]*((y-y0[3])**2))
        dV_dy = A[0]*(c[0]*2*(y-y0[0]) + b[0]*(x-x0[0]))*exp(a[0]*((x-x0[0])**2) + b[0]*(x-x0[0])*(y-y0[0]) + c[0]*((y-y0[0])**2)) + A[1]*(c[1]*2*(y-y0[1]) + b[1]*(x-x0[1]))*exp(a[1]*((x-x0[1])**2) + b[1]*(x-x0[1])*(y-y0[1]) + c[1]*((y-y0[1])**2)) + A[2]*(c[2]*2*(y-y0[2]) + b[2]*(x-x0[2]))*exp(a[2]*((x-x0[2])**2) + b[2]*(x-x0[2])*(y-y0[2]) + c[2]*((y-y0[2])**2)) + A[3]*(c[3]*2*(y-y0[3]) + b[3]*(x-x0[3]))*exp(a[3]*((x-x0[3])**2) + b[3]*(x-x0[3])*(y-y0[3]) + c[3]*((y-y0[3])**2)) + A[3]*(a[3]*2*(x-x0[3]) + b[3]*(y-y0[3]))*exp(a[3]*((x-x0[3])**2) + b[3]*(x-x0[3])*(y-y0[3]) + c[3]*((y-y0[3])**2))
    
    # Force calculation
    fx = -dV_dx - gamma*px + np.random.normal(0,std_dev)
    fy = -dV_dy - gamma*py + np.random.normal(0,std_dev)
    return fx,fy
 
def vv_step(phasepoint,dt,beta,gamma,pes_type):
    x = phasepoint[0]
    y = phasepoint[1]
    px = phasepoint[2]
    py = phasepoint[3]
    fx = phasepoint[4]
    fy = phasepoint[5]
    px = px + (1/2)*dt*fx
    py = py + (1/2)*dt*fy
    x = x + dt*px
    y = y + dt*py
    fx,fy = force(x,y,px,py,dt,beta,gamma,pes_type)
    px = px + (1/2)*dt*fx
    py = py + (1/2)*dt*fy
    return np.asarray([x,y,px,py,fx,fy])
 
def calc_op(op_type,x,y):
    if op_type == 1:
        return x
    elif op_type == 2:
        return y
    elif op_type == 3:
        return x + y
    else:
        sys.exit("Invalid choice of OP")

