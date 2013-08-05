'''
Created on Jul 8, 2013

@author: lee
'''
from mean_shift import *
from numpy import linalg

numberofdataPoints = 1000
numberofclusers = 12
dimensions = 2000
k=100


if __name__ == '__main__':
    pass
    
    d = dimensions
    pts,cntrs =  getDataPoints(numberofdataPoints,dimensions,numberofclusers)
    
    prs = random.randn(d,2)#/array([[1./2. , 1./2. ] for i in range(d)])
    
    
    
    prsi = linalg.inv
    
    K = dot(rPts,prs)
    
    V = doMeanShift(K)