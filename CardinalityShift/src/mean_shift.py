import random
from pyflann import *
#import random

k=10
def getDataPoints(part,d,clu):
    ret = [[] for j in xrange(part*clu)]
    clusterCenters = []
    for i in range(clu):
        variance = .5+random.random()*4
        means = [random.random()*10 for b in range(d)] 
        clusterCenters.append(means)
        for j in range(part):
            p = [0]*d
            for k in xrange(d):
                pts = random.random()
                p[k]= (pts*variance+means[k])
            ret[i*part+j]=p
    return ret,clusterCenters
    
def dist(a,b):return sum((a[i]-b[i])**2 for i in xrange(len(a)))

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

def getBandwidth(vec, d):
    '''
    return distance of kth nearest neighbor
    '''
    l = sorted(d)
    return l[k-1]

def prune(vec,d):
    '''
        this is where lsh db-gen would be
    '''
    from pylab import array
    neighbors,nearest = flann.nn(array(d),array(vec),k)

    ret = []
    for idx in neighbors[0]:
        ret.append(d[idx])
    
    
    '''
    l =sqdistance(vec,d)
    dist = sorted( [ (l[i] , i)   for i in xrange(len(d))] ) 
    ret = [[]]*k
    for i in xrange(k):
        ret[i] = d[ dist[i][1] ]

    '''
    
    
    return ret
    
import copy
from math import exp
flann = FLANN()
def doMeanShift(dataPoints):
    '''
        perform the mean shift operations
    '''
    from pylab import array
    origData = copy.deepcopy(dataPoints)
    
    if k != 0: 
        
        flann.build_index(array(dataPoints))
    
    threshold = 1e-6
    newVecs = []

    for vec in dataPoints:
        diff = threshold+1

        if k!=0:
            pruned = prune(vec,origData)
        else:
            pruned =origData# prune(vec,origData)       
        
        while diff > threshold:
            d = sqdistance(vec,pruned)
            bandwidth= getBandwidth(vec,d)
            kernelDist = [exp(-dd)/(bandwidth*bandwidth) for dd in d] 
            numerator = [0.0]*len(vec)      
            for i in xrange(len(vec)):
                for j in xrange(len(kernelDist)):
                    numerator[i] +=  kernelDist[j] * origData[j][i]      
                       
            denom = sum(kernelDist)          
            tmp = [num/denom for num in numerator] 
            diff = max([abs(tmp[i]-vec[i]) for i in xrange(len(vec)) ] )
            vec = tmp
        
        #only add unique vectors
        ct = 0
        while ct< len( newVecs) and dist(vec,newVecs[ct])>threshold :
            ct+=1
            
        if ct== len(newVecs):         
                newVecs.append(vec)   
                
    return newVecs
    
def zero(m,n):
    # Create zero matrix
    new_matrix = [[0 for row in xrange(n)] for col in xrange(m)]
    return new_matrix

def mult(matrix1,matrix2):
    # Matrix multiplication
    if len(matrix1[0]) != len(matrix2):
        # Check matrix dimensions
        print 'Matrices must be m*n and n*p to multiply!'
    else:
        # Multiply if correct dimensions
        new_matrix = zero(len(matrix1),len(matrix2[0]))
        for i in xrange(len(matrix1)):
            for j in xrange(len(matrix2[0])):
                for k in xrange(len(matrix2)):
                    new_matrix[i][j] += matrix1[i][k]*matrix2[k][j]
        return new_matrix
 


def drawPts(V,pts):
    
    
    x = [p[0] for p in V]
    y = [p[1] for p in V]
    z = [p[2] for p in V]
    
    

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    
    fig = plt.figure(1)
    ax = Axes3D(fig)
    
    colorpalette = ['blue','green','red','purple','yellow','orange','gray','brown','magenta','cyan','white']
    clucolors = ['']*len(pts)
    clu = [V[0]]
    for p in range(len(V)):
        lst = 10000.  
        for c in range(len(clu)):
            if dist(V[p],clu[c]) <lst:
                lst = dist(V[p],clu[c])
                clucolors[p]=colorpalette[c % len(colorpalette)]    
        if lst > 1e-2:
            clucolors[p]=colorpalette[ len(clu) % (len(colorpalette))]
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
                clucolors[p]=colorpalette[c% len(colorpalette)]
    ax.scatter(x,y,z, color=clucolors,marker='o')
    


def randProjectPts(pts,cluDim):
    d = len(pts[0])
    r = 1.0/float(len(pts[0])**.5)
    M= ([[random.normalvariate(0,1)*r for a in range(cluDim)] for b in range(d)])

    return mult(pts,M)
    
if __name__ == '__main__':
    pass        

    from random import shuffle
    
    numberofdataPoints = 1000
    numberofclusers = 5
    dimensions = 3
        
    pts,cntrs =  getDataPoints(numberofdataPoints,dimensions,numberofclusers)
    #shuffle(pts) //not really needed to prove algorithm correctness
    
    
    print cntrs
    V = doMeanShift(pts)
    #print V
     
    print "projecting data points\n"
    #ptsRedux = randProjectPts(pts,3)
    #print len(ptsRedux),len(ptsRedux[0])
    print "performing mean shift on projected data\n"
    
    #V = doMeanShift(ptsRedux)
    print V
    
    

