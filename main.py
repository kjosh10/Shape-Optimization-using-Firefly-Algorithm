"""
importing necessary dependencies
"""

import numpy as np
import matplotlib.pyplot as plt
import sys,os
sys.path.append(os.getcwd()+"/DO_code")
from util.CST_FEM import *
from util.methods import *
from util.FireFlyAlgorithm import *
from util.MeshDomain import *

#%%

lx = 5.5
ly = 0.2
h  = 0.1
domain = [(0,0),(lx,ly)]
n_x = 21
n_y = 5
E = 2.1e11
nu = 0.2
  
parameters = [5, 10, 1, 0.2, 1]
# [no.FireFlies, maxIteration, alpha, betaMin, gamma]


#%%
fireFly = FireFly(parameters,domain,n_x,n_y,E,nu,h)
steps = fireFly.runAlgorithm()
min_index = np.argmin(fireFly.volume_)
bestStep = steps[min_index]
fireFly.stress_[min_index]
fireFly.bestStep
fireFly.bestStress
fireFly.bestVolume

#np.arange(fireFly.maxIteration)
fig, ax = plt.subplots()
ax.plot(np.arange(1,fireFly.maxIteration+1), fireFly.bestStep)
ax.set_xlabel("Iteration")
ax.set_ylabel("Curve height [m]")
ax.set_title("Curve height change with iterations")
fig.show()
plt.savefig("CurveHeight.png",dpi=800)

fig, ax = plt.subplots()
ax.plot(np.arange(1,fireFly.maxIteration+1), fireFly.bestStress)
ax.set_xlabel("Iteration")
ax.set_ylabel("Stress [MPa]")
ax.set_title("Stress change with iterations")
fig.show()
plt.savefig("Stress.png",dpi=800)

fig, ax = plt.subplots()
ax.plot(np.arange(1,fireFly.maxIteration+1), fireFly.bestVolume)
ax.set_xlabel("Iteration")
ax.set_ylabel("Volume [m^3]")
ax.set_title("Volume change with iterations")
fig.show()
plt.savefig("Volume.png",dpi=800)

fig, ax = plt.subplots()
ax.plot(np.arange(1,fireFly.no_fireflies+1), fireFly.stepInitialized, label="Initialized Curve heigths")
ax.plot(np.arange(1,fireFly.no_fireflies+1), steps, label="Final Curve heigths")
ax.set_xlabel("No. of Fireflies")
ax.set_ylabel("Curve Height")
ax.set_title("Curve height change with iterations after random Initialization")
ax.legend()
fig.show()
plt.savefig("CurveHeight Initialized.png",dpi=800)

#%%


volume_ = []
stress_ = []
steps_  = []
k = 1
for step  in np.array([bestStep]):
    
    y_mid = step
    points = [(0,0),(lx,0)]
    curve_y_mid = bezierCurveTopPoint(y_mid,points)
    points_array_Curve = [(0,0),(lx/2,curve_y_mid),(lx,0)]
    meshDomain_Beam = MeshDomain(domain, points_array_Curve, n_x, n_y)
    


    #number of nbc

    Xs, Ys, pts = MeshDomain_Beam.meshPointsCalc()
    
    cst_arr = np.empty(len(Xs),dtype=object)
    
    
    #distributed load
    n_nbc = n_x
    nbc_arr = np.empty(n_nbc,dtype=object)
    R = -12000000
    nodal_force = R/(n_x)
    for i in range(n_x):
        nbc_arr[i] = Nbc(x=i*lx/(n_x-1),y=ly,fx=0,fy=R/(n_x))
        
    #element array
    for i in range(len(Xs)):
        cst_arr[i] = CST(Xs[i],Ys[i],E,nu)
        cst_arr[i].h = h
        
        
        
    #dbc array
    points_grid_Mesh = MeshDomain_Beam.points_grid
    points_grid_Mesh = np.array(points_grid_Mesh) 
    points_grid_Mesh.shape
    dbc_Points = []
    for i in range(points_grid_Mesh.shape[0]):
        dbc_Points.append(points_grid_Mesh[i][0])
        dbc_Points.append(points_grid_Mesh[i][-1])
    dbc_Points.append([0,0])
    dbc_Points.append([lx,0])
    
    dbc_arr = np.empty(len(dbc_Points),dtype=object)
    for i in range(len(dbc_Points)):
        dbc_arr[i] = Dbc(dbc_Points[i][0],dbc_Points[i][1],1,1)
    
    
    model1 = Model(cst_arr,dbc_arr,nbc_arr)
    model1.solve()
    model1.plot()
    cst_arr = model1.update_elements()
    
    print("-------STEP------", step)
    print("total volume = ", model1.volume())
    print("max stress = " , model1.max_stress())
    print("min disp = ", model1.min_u())
    u = model1.u
    stress_.append(model1.max_stress()/1e6)
    volume_.append(model1.volume())
    steps_.append(k)

    k += 1


#%%

lx = 5
ly = 0.2
h  = 0.1
domain = [(0,0),(lx,ly)]
n_x = 50
n_y = 5
E = 2.1e11
nu = 0.2

meshDomain_Beam = MeshDomain(domain, n_x=n_x, n_y=n_y)
    


 #number of nbc

Xs, Ys, pts = MeshDomain_Beam.meshPointsCalc()
    
cst_arr = np.empty(len(Xs),dtype=object)
    
    
#distributed load
n_nbc = n_x
nbc_arr = np.empty(n_nbc,dtype=object)
R = -1000000
# R = -12000
nodal_force = R/(n_x)
for i in range(n_x):
    nbc_arr[i] = Nbc(x=i*lx/(n_x-1),y=ly,fx=0,fy=R/(n_x))
        
 #element array
for i in range(len(Xs)):
    cst_arr[i] = CST(Xs[i],Ys[i],E,nu)
    cst_arr[i].h = h
        
        
        
#dbc array
points_grid_Mesh = MeshDomain_Beam.points_grid
points_grid_Mesh = np.array(points_grid_Mesh) 
points_grid_Mesh.shape
dbc_Points = []
for i in range(points_grid_Mesh.shape[0]):
    dbc_Points.append(points_grid_Mesh[i][0])
    dbc_Points.append(points_grid_Mesh[i][-1])
# dbc_Points.append([0,0])
# dbc_Points.append([lx,0])
    
dbc_arr = np.empty(len(dbc_Points),dtype=object)
for i in range(len(dbc_Points)):
    dbc_arr[i] = Dbc(dbc_Points[i][0],dbc_Points[i][1],1,1)

model1 = Model(cst_arr,dbc_arr,nbc_arr)
model1.solve()
model1.plot()
cst_arr = model1.update_elements()
print("total volume = ", model1.volume())
print("max stress = " , model1.max_stress())
print("min disp = ", model1.min_u())

#%%  //plot1

plt.plot(steps_,stress_,"ro",ms=2)
plt.xlabel("step")
plt.ylabel("stress (mPa)")

#%%

plt.plot(steps_,volume_,"bo",ms=2)
plt.xlabel("step")
plt.ylabel("volume (m3)")




