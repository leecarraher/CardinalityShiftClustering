from numpy import random
numberofdataPoints = 100
numberofclusers = 3
dimensions = 3
def getDataPoints(part,d,clu):
    ret = [[] for j in xrange(part*clu)]
    for i in range(clu):
        variance = random.random()*4
        means = [random.random(1)*10 for b in range(d)] 
        
        for j in range(part):
            p = [0]*d
            for k in xrange(d):
                pts = random.random()
                p[k]= (pts*variance+means[k])[0]
            ret[i*part+j]=p
    return ret

pts =  getDataPoints(numberofdataPoints,dimensions,numberofclusers)
   
    
def sqdistance(a,b):
    '''
        a is a vector
        b is a matrix
    '''
    aa = 0
    for i in xrange(len(a)):
        aa = aa + a[i] * a[i]
    
    bb = [aa]*len(b)
    
    for i in xrange(len(b)):
        
    
        for j in xrange(len(b[i])):
            bb[i] = bb[i]+b[i][j]*b[i][j]
        for j in xrange(len(a)):
            bb[i] = bb[i] - 2*(a[j] * b[i][j])
    
    return bb
    
    
def getBandwidth(curPt):
    '''
    return distance of kth nearest neighbor




    function bandwidth = getBandWith(curDataPoint,origDataPoints,euclideanDist,useKNNToGetH)
            if (useKNNToGetH == false)
                    bandwidth = H;
            else:
                    [sortedEuclideanDist,indexInOrig] = sort(euclideanDist);
                    bandwidth = sortedEuclideanDist(k);
            end
    end
    
    '''
    return 1.0
#print sqdistance([1,2,3],[[1,2,4],[3,1,2],[4,1,2],[4,1,2],[4,1,2],[4,5,7]])

import copy
from math import exp
def doMeanShift(dataPoints):
    '''
        perform the mean shift operations
    '''

    origData = copy.deepcopy(dataPoints)
    
    threshold = 1e-3
    newVecs = []
    for vec in dataPoints:
        diff = threshold+1
        while diff > threshold:
            d = sqdistance(vec,origData)
            bandwidth= getBandwidth(vec)
            kernelDist = [exp(-dd)/(bandwidth*bandwidth) for dd in d]
            
            numerator = [0.0]*len(vec)
            
            for i in xrange(len(vec)):
                for j in xrange(len(kernelDist)):
                    numerator[i] +=  kernelDist[j] * origData[j][i]
            
            denom = sum(kernelDist)
            
            tmp = [num/denom for num in numerator]
            
            diff = max([abs(tmp[i]-vec[i]) for i in xrange(len(vec)) ] )
            vec = tmp
        
        newVecs.append(vec)
        
        print diff
    return newVecs
        
V = doMeanShift(pts)

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


fig = plt.figure(1)
ax = Axes3D(fig)

x = [p[0] for p in pts]
y = [p[1] for p in pts]
z = [p[2] for p in pts]

ax.scatter(x,y,z, marker='o')

x = [p[0] for p in V]
y = [p[1] for p in V]
z = [p[2] for p in V]

ax.scatter(x,y,z, marker='x')
plt.show()





