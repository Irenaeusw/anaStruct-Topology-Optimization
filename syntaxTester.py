import anaTool 
import numpy as np 
from anastruct import SystemElements 
import timeit 
import matplotlib.pyplot as plt 

xHexagons = 9
yHexagons = 17
centerPoint = (0,0) 

saveDir = r"C:\Users\Irena\Desktop\Zurob Research\honeycomb lattice anastruct"

nodeDistance = 2.5 # In millimeters
latticeParam = np.round( (nodeDistance/2)/np.tan(np.pi/6) , 4)  
force = -1000 # In newtons 
elasticModulus = 17000 # In pascals 
moment = elasticModulus/3 
area = 0.5 # Substitute for wall thickness, in mm 

ss = anaTool.initializeLatticeShell( startPoint=(0,0), latticeParam=latticeParam, xHexagons=xHexagons, yHexagons=yHexagons, 
        force=force, theta=0,  elasticModulus=elasticModulus, moment=moment, area=area, saveDir=saveDir)

'''
E - Young's Modulus
A = Area
I = Moment of Inertia 
q
EA - Standard axial stiffness of element 
EI - Standard bending stiffness of elements 
'''







