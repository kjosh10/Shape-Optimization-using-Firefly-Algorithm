"""
importing necessary dependencies
"""

import numpy as np
import matplotlib.pyplot as plt

class CST:
    """
    this script is used to create a CST (Constant Strain Triangle) Element

    Attributes:
        x: list, a list of all the x co-ordinates
        y: list, a list of all the y co-ordinates 
        E: float, young's modulus of the material
        nu: float, poisson's ratio of the material
    """
    def __init__(self,x,y,E,nu):
        self.x = [round(num, 3) for num in x]
        self.y = [round(num, 3) for num in y]
        self.E = E
        self.nu = nu
        self.h = 1.0
        
        self.u = np.zeros(6)
        
        
        self.strain = np.zeros(3)
        
        
        self.stress = np.zeros(3)
       
        #area (A)
        self.A = (1.0/2.0)*x[1]*y[2] - 1.0/2.0*x[1]*y[3] - 1.0/2.0*x[2]*y[1] + (1.0/2.0)*x[2]*y[3] + (1.0/2.0)*x[3]*y[1] - 1.0/2.0*x[3]*y[2]
        self.A = abs(self.A)
        
        #derivative of shape functions (B)
        self.B = np.array([[         y[2] - y[3],                0,     -y[1] + y[3],                0,  y[1] - y[2],           0],
                                [0, -x[2] + x[3],                0,      x[1] - x[3],            0, -x[1] + x[2]],
                     [-x[2] + x[3],  y[2] - y[3],      x[1] - x[3],     -y[1] + y[3], -x[1] + x[2],  y[1] - y[2]]])
        self.B = self.B * 1/(2*self.A)
        #Constitutive matrix (C)
        self.C = np.zeros(9).reshape(3,3)
        self.C[0][0] = 1; self.C[0][0] = 1; self.C[0][1]=nu; self.C[1][0] = nu; self.C[1][1]=1;  self.C[2][2] = 0.5*(1-nu);
        
        self.C *= E/(1-nu*nu)
        
        #element stiffness matrix (K)
        self.K = self.B.transpose()@self.C@self.B
        self.K = self.h*self.A*self.K
        
    def plot(self,scale=0):
        """
        method to plot the CST Element

        :param scale: float, float to create the scaling while plotting
        """
        plt.plot(self.x[1:]+scale*self.u[0:3],self.y[1:]+scale*self.u[3:6],"b-",linewidth=0.3)
        plt.plot((self.x[1]+scale*self.u[0],self.x[3]+scale*self.u[2]),(self.y[1]+scale*self.u[3],self.y[3]+scale*self.u[5]),"b-",linewidth=0.3)
        
        
        
    #calculate stress
    def calc_stress(self):
        """
        method to compute stress in the CST Element

        :return stress: float, stress in the element
        """
        #strain
        self.eps = self.B @ self.u
        #stress
        self.sigma = self.C @ self.eps
        sx = self.sigma[0]
        sy = self.sigma[1]
        txy = self.sigma[2]
        stress = np.sqrt(sx**2-(sx*sy)+sy**2+3*txy**2)
        return stress


class Dbc:
    """
    this script is used to create a Direchlet Boundary condition (DBC)

    Attributes:
        x: float, x co-ordinate
        y: float, y co-ordinate
        ux: boolean, a boolean to say the presence of DBC in x direction
        uy: boolean, a boolean to say the presence of DBC in y direction
    """
    def __init__(self,x,y,ux,uy):
        self.x = round(x,3)
        self.y = round(y,3)
        self.ux = ux     #boolean
        self.uy = uy     #boolean

class Nbc:
    """
    this script is used to create a Neumann Boundary condition (NBC)

    Attributes:
        x: float, x co-ordinate
        y: float, y co-ordinate
        fx: float, force value in x direction
        fy: boolean, force value in y direction
    """
    def __init__(self,x,y,fx,fy):
        self.x = round(x,3)
        self.y = round(y,3)
        self.fx = fx
        self.fy = fy

class Model:
    """
    this script is used to create a Neumann Boundary condition (NBC)

    Attributes:
        cst_arr: list, a list of all the cst elements
        dbc_arr: list, a list of all the DBC
        nbc_arr: list, a list of all the NBC
    """
    def __init__(self,cst_arr,dbc_arr,nbc_arr):
        self.cst_arr=cst_arr
        self.dbc_arr=dbc_arr
        self.nbc_arr=nbc_arr
        
        #create unique point array
        unique_node_set = set()
        coordx = 0; coordy=0;
        for i in range(len(cst_arr)):
            coordx=round(self.cst_arr[i].x[1],3); coordy=round(self.cst_arr[i].y[1],3); unique_node_set.add((coordx,coordy));
            coordx=round(self.cst_arr[i].x[2],3); coordy=round(self.cst_arr[i].y[2],3); unique_node_set.add((coordx,coordy));
            coordx=round(self.cst_arr[i].x[3],3); coordy=round(self.cst_arr[i].y[3],3); unique_node_set.add((coordx,coordy));
    
        self.unique_node_list = list(unique_node_set)
        
        #print(self.unique_node_list)
        
        self.Kglob = np.zeros((2*len(self.unique_node_list),2*len(self.unique_node_list)))
        self.u = np.zeros(2*len(self.unique_node_list))
        self.f = np.zeros(2*len(self.unique_node_list))
        #assembly to global stiffness matrix
        #find the index of point in set 
        for i in range(self.cst_arr.size):
            point1 = (self.cst_arr[i].x[1],self.cst_arr[i].y[1])     #(x,y)
            point2 = (self.cst_arr[i].x[2],self.cst_arr[i].y[2])
            point3 = (self.cst_arr[i].x[3],self.cst_arr[i].y[3])
                
            p1_idx = self.unique_node_list.index(point1)
            p2_idx = self.unique_node_list.index(point2)
            p3_idx = self.unique_node_list.index(point3)
            
            kelement = cst_arr[i].K
            k11 = kelement[0:2,0:2]; k12 = kelement[0:2,2:4]; k13 = kelement[0:2,4:6];
            k21 = kelement[2:4,0:2]; k22 = kelement[2:4,2:4]; k23 = kelement[2:4,4:6];
            k31 = kelement[4:6,0:2]; k32 = kelement[4:6,2:4]; k33 = kelement[4:6,4:6];
            
            self.Kglob[2*p1_idx:2*p1_idx+2,2*p1_idx:2*p1_idx+2] += k11
            self.Kglob[2*p1_idx:2*p1_idx+2,2*p2_idx:2*p2_idx+2] += k12
            self.Kglob[2*p1_idx:2*p1_idx+2,2*p3_idx:2*p3_idx+2] += k13
            
            self.Kglob[2*p2_idx:2*p2_idx+2,2*p1_idx:2*p1_idx+2] += k21
            self.Kglob[2*p2_idx:2*p2_idx+2,2*p2_idx:2*p2_idx+2] += k22
            self.Kglob[2*p2_idx:2*p2_idx+2,2*p3_idx:2*p3_idx+2] += k23
            
            self.Kglob[2*p3_idx:2*p3_idx+2,2*p1_idx:2*p1_idx+2] += k31
            self.Kglob[2*p3_idx:2*p3_idx+2,2*p2_idx:2*p2_idx+2] += k32
            self.Kglob[2*p3_idx:2*p3_idx+2,2*p3_idx:2*p3_idx+2] += k33
        
        
        #force vector
        for i in range(len(nbc_arr)):
            point = (nbc_arr[i].x,nbc_arr[i].y)
            idx = self.unique_node_list.index(point)
            self.f[2*idx]=nbc_arr[i].fx
            self.f[2*idx+1]=nbc_arr[i].fy
            
        #impose dbc
        for i in range(len(dbc_arr)):
            point = (dbc_arr[i].x,dbc_arr[i].y)
            idx = self.unique_node_list.index(point)
            if(dbc_arr[i].ux==1):
                self.Kglob[2*idx][2*idx]*=10e12
            if(dbc_arr[i].uy==1):
                self.Kglob[2*idx+1][2*idx+1]*=10e18
    
    def solve(self):
        """
        method to solve the model using FEM method

        :param u: np.ndarray, an array of the solution of displacements
        """
        self.u = np.linalg.solve(self.Kglob,self.f)
        
    def update_elements(self):
        """
        method to update all the elements with the solution

        :param cst_arr: list, a list of updated cst elements
        """
        for i in range(len(self.cst_arr)):
            for j in range(1,4):
                idx = self.unique_node_list.index((self.cst_arr[i].x[j],self.cst_arr[i].y[j]))
                self.cst_arr[i].u[2*j-2] = self.u[2*idx]   
                self.cst_arr[i].u[2*j-1] = self.u[2*idx+1]
            
        return self.cst_arr                    

    
    def volume(self):
        """
        method to compute volume of the model

        :param v: float, volume of the model
        """
        v = 0.0;
        for i in range(len(self.cst_arr)):
            v += self.cst_arr[i].A * self.cst_arr[i].h
        return v


    def plot(self,scale=0):
        """
        method to plot the model
        """
        for i in range(len(self.cst_arr)):
            self.cst_arr[i].plot(scale)
            plt.axis("equal")
        
    def max_stress(self):
        """
        method to compute maximum stress in the model

        :return max_stress: float, maximum stress in the model 
        """
        max_stress = 0
        for i in range(len(self.cst_arr)):
            if(abs(self.cst_arr[i].calc_stress())>max_stress):
                max_stress = abs(self.cst_arr[i].calc_stress())
        return max_stress
    
    
    def min_u(self):
        """
        method to compute minimum displacement 

        :return min_u: float, minimum displacement 
        """
        min_u = 0.
        for i in range(len(self.cst_arr)):
            u_ = self.cst_arr[i].u
            if(u_[1]<min_u):
                min_u = u_[1]
            if(u_[3]<min_u):
                min_u = u_[3]
            if(u_[5]<min_u):
                min_u = u_[5]
        return min_u
#%%
#axial loaded plane-stress
if(__name__=="__main__"):

    E = 2.1e11
    nu = 0.0
    cst_arr = np.empty(16,dtype=object)
    dbc_arr = np.empty(2,dtype=object)
    nbc_arr = np.empty(2,dtype=object)
    #serie 1 elements
    for i in range(8):
        cst_arr[i] = CST(x=(0,4*i,4*i+4,4*i+4),
                         y=(0,0,0,3),       
                         E=E,nu=nu)
    #series 2 elements
    for i in range(8):
        cst_arr[i+8] = CST(x=(0,4*i,4*i,4*i+4),
                           y=(0,0,3,3),
                           E=E,nu=nu)
    
    #DBC
    dbc_arr[0] = Dbc(0,0,1,1)
    dbc_arr[1] = Dbc(0,3,1,1)
    
    #NBC
    nbc_arr[0] = Nbc(32,0,500000,0)
    nbc_arr[1] = Nbc(32,3,500000,0)    
    
    #solve
    model1 = Model(cst_arr,dbc_arr,nbc_arr)
    model1.solve()
    cst_arr = model1.update_elements()
    cst_arr[1].calc_stress()
    
    model1.u
    model1.plot()
    model1.plot(scale=10000)
    
    model1.volume()
    
    #analytical
    print(500000*32/E/3)  #FL/EA
