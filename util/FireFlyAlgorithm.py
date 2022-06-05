"""
importing necessary dependencies
"""

import numpy as np
import random as rd
import math
from CST_FEM import *
from DO import *
import sys

#%%

class FireFly:
    """
    this script is used to implement Firefly algorithm

    Attributes:
        parameters: list, a list of all the parameters
        n_x: int, total co-ordinates in x direction
        n_y: int, total co-ordinates in y direction
        E: float, young's modulus of the material
        nu: float, poisson's ratio of the material
        h: float, thickness of the model
    """
    def __init__(self,parameters,domain,n_x,n_y,E,nu,h): 
        self.no_fireflies = parameters[0]
        self.maxIteration = parameters[1]
        self.alpha = parameters[2]
        self.betaMin = parameters[3]
        self.gamma = parameters[4]
        self.domain = domain
        self.n_x = n_x
        self.n_y = n_y
        self.E = E
        self.nu = nu
        self.h = h
        self.lowerBound = domain[0][1]
        self.upperBound = domain[1][1] -1e-3   # To avoid sharp corners
        
        # Total number of function evaluations
        self.numEvaluation = self.no_fireflies*self.maxIteration
        # Dimension 
        self.dim = 1   # Only one design variable
        
        # Initial values of an array  ????
        #self.zn=ones(n,1)*10^100;
        
        
        
    #-----------The initial locations of n fireflies------------#
    def initialize_Fireflies(self):
        """
        method to initialize the fireflies

        :return steps: list, a list of all the firefly population
        """
        steps = np.zeros(self.no_fireflies)
        for i in range(self.no_fireflies):
            steps[i] = self.lowerBound + rd.random()*(self.upperBound - self.lowerBound)
        
        return steps
                             
    def eval(self, steps):
        """
        method to evaluate the performance of the fireflies

        :param steps: list, a list of all the firefly population
        :return volume_, float, volume of the model
        """
        self.volume_ = []
        self.stress_ = []
        self.steps_  = []
        k = 1
        lx = self.domain[1][0]
        ly = self.domain[1][1]
        
        for step in steps:
            
            y_mid = step
            points = [(0,0),(lx,0)]
            y_mid_Bezier = bezierCurveTopPoint(y_mid,points)
            points_array_Curve = [(0,0),(lx/2,y_mid_Bezier),(lx,0)]
            meshDomain_Beam = meshDomain(self.domain, points_array_Curve, self.n_x, self.n_y)
            
            #number of nbc
        
            Xs, Ys, pts = meshDomain_Beam.meshPointsCalc()
            
            cst_arr = np.empty(len(Xs),dtype=object)
            
            
            #distributed load
            n_nbc = self.n_x
            nbc_arr = np.empty(n_nbc,dtype=object)
            R = -12000000
            nodal_force = R/(self.n_x)
            for i in range(self.n_x):
                nbc_arr[i] = Nbc(x=i*lx/(self.n_x-1),y=ly,fx=0,fy=R/(self.n_x))
                
            #element array
            for i in range(len(Xs)):
                cst_arr[i] = CST(Xs[i],Ys[i],self.E,self.nu)
                cst_arr[i].h = self.h
                
                
                
            #dbc array
            points_grid_Mesh = meshDomain_Beam.points_grid
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
            #model1.plot()
            cst_arr = model1.update_elements()
            
            self.stress_.append(model1.max_stress()/1e6)
            self.volume_.append(model1.volume())
        
        #-------------- Constraint ---------------------------------#
        yieldStress = 350 # MPa
               
        for i in range(self.no_fireflies):
           penalties = self.penaltyCheckerNApplier(np.array(self.stress_), yieldStress)
        
        #-------------- Computing Objective Function ---------------#
        self.volume_ = np.array(self.volume_) + penalties
        
        if (steps.shape[0] == 1):
            return self.volume_
    
    
    
    def penaltyCheckerNApplier(self, stresses, yieldStress):
        """
        method to check and apply the penalty for the optimization problem

        :param stresses: list, stresses of the model
        :param yieldStress: float, yield stress of the material
        """
        # -------------Define Penalty --------------------#
        penalty = 1e10
        
        stresses = np.array(stresses)
        yieldStress = np.ones(stresses.shape[0])*yieldStress
        g = stresses - yieldStress
        penalties = np.zeros(g.shape[0])
        for index, g_val in enumerate(g):
            if (g_val > 0):
                penalties[index] = penalty*g_val
        
        return penalties
        
        
    def runAlgorithm(self):
        """
        method to run the Firefly algorithm

        :return steps: list, a list of all the firefly population
        """
        # Creation of Population
        steps = self.initialize_Fireflies()
        self.stepInitialized = np.copy(steps)
        self.bestStress = np.zeros(self.maxIteration)
        self.bestVolume = np.zeros(self.maxIteration)
        self.bestStep = np.zeros(self.maxIteration)
        # print(steps)
        #self.eval(steps)
        self.scale = self.upperBound - self.lowerBound
        print("Loop Started with one dot equal to 5 Iterations:")
        for iter in range(self.maxIteration):
            for i in range(self.no_fireflies):
                for j in range(self.no_fireflies):
                    if (self.eval(np.array([steps[i]])) <= self.eval(np.array([steps[j]]))):
                        steps[i] = steps[i]
                    else:
                        r = math.sqrt((steps[j]-steps[i])**2)
                        beta = self.betaMin*math.exp(-self.gamma*r**2)
                        increment_step = self.alpha*(rd.random()-0.5)*self.scale
                        stepNew = steps[i] + beta*(steps[j]-steps[i])+increment_step
        
                        # Check Bounds
                        if (stepNew > self.upperBound):
                            stepNew = self.upperBound
                        elif (stepNew < self.lowerBound):
                                stepNew = self.lowerBound
                    
                        # Greedy selection
                        if (self.eval(np.array([stepNew])) < self.eval(np.array([steps[i]]))):
                            steps[i] = stepNew;
            self.eval(steps)
            min_index = np.argmin(self.volume_)
            self.bestVolume[iter] = self.volume_[min_index]
            self.bestStep[iter] = steps[min_index]
            self.bestStress[iter] = self.stress_[min_index]
            if (((iter+1) % 5 == 0) and iter != 0):
                sys.stdout.write(".")
                sys.stdout.flush()
                
        self.eval(steps)
        
        return steps
               
    
