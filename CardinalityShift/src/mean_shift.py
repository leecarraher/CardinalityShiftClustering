from numpy import random
numberofdataPoints = 50
numberofclusers = 5
dimensions = 3
k=50
def getDataPoints(part,d,clu):
    ret = [[] for j in xrange(part*clu)]
    clusterCenters = []
    for i in range(clu):
        variance = .3+random.random()*5
        means = [random.random(1)*10 for b in range(d)] 
        clusterCenters.append(means)
        for j in range(part):
            p = [0]*d
            for k in xrange(d):
                pts = random.random()
                p[k]= (pts*variance+means[k])[0]
            ret[i*part+j]=p
    return ret,clusterCenters
    
def sqdistance(a,b):
    '''
        a is a vector
        b is a matrix
        this is basically a distances list, and should be pruned by lsh-knn
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

def getBandwidth(vec, d,k):
    '''
    return distance of kth nearest neighbor
    '''
    return d[k]

import copy
from math import exp
def doMeanShift(dataPoints):
    '''
        perform the mean shift operations
    '''
    origData = copy.deepcopy(dataPoints)
    
    threshold = 1e-5
    newVecs = []
    for vec in dataPoints:
        diff = threshold+1
        while diff > threshold:
            d = sqdistance(vec,origData)
            
            #starts
            d = sorted(d)
            bandwidth= getBandwidth(vec,d,int(k))
            
            kernelDist = [exp(-dd)/(bandwidth*bandwidth) for dd in d[:k]]
            numerator = [0.0]*len(vec)
            
            for i in xrange(len(vec)):
                for j in xrange(len(kernelDist)):
                    numerator[i] +=  kernelDist[j] * origData[j][i]
            '''
            kernelDist = [exp(-dd)/(bandwidth*bandwidth) for dd in d] 
            numerator = [0.0]*len(vec)      
            for i in xrange(len(vec)):
                for j in xrange(len(kernelDist)):
                    numerator[i] +=  kernelDist[j] * origData[j][i]         
            '''
            denom = sum(kernelDist)          
            tmp = [num/denom for num in numerator] 
            diff = max([abs(tmp[i]-vec[i]) for i in xrange(len(vec)) ] )
            vec = tmp
        newVecs.append(vec)   
        print diff
    return newVecs
        
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure(1)
ax = Axes3D(fig)

pts,cntrs =  getDataPoints(numberofdataPoints,dimensions,numberofclusers)

V = doMeanShift(pts)
x = [p[0] for p in V]
y = [p[1] for p in V]
z = [p[2] for p in V]

def dist(a,b):return sum((a[i]-b[i])**2 for i in xrange(len(a)))

colorpalette = ['blue','green','red','purple','yellow','orange','gray','brown','magenta','cyan','white']
clucolors = ['']*len(pts)
clu = [V[0]]
for p in range(len(V)):
    lst = 10000.  
    for c in range(len(clu)):
        if dist(V[p],clu[c]) <lst:
            lst = dist(V[p],clu[c])
            clucolors[p]=colorpalette[c]    
    if lst > 1e-3:
        clucolors[p]=colorpalette[len(clu)]
        clu.append(V[p])
ax.scatter(x,y,z, color=clucolors, marker='x')
x = [p[0] for p in pts]
y = [p[1] for p in pts]
z = [p[2] for p in pts]

clucolors = ['']*len(pts)
for p in range(len(V)):
    lst = 10000.
    
    for c in range(len(clu)):
        if dist(V[p],clu[c]) <lst:
            lst = dist(V[p],clu[c])
            clucolors[p]=colorpalette[c]
ax.scatter(x,y,z, color=clucolors,marker='o')

print "Original Centers"
cntrs = [ [ cc[0] for cc in c]   for c in cntrs]
print cntrs
print "Predicted Centers"
print clu


plt.show()


