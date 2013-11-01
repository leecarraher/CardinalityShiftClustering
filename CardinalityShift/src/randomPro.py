
from pylab import *
   
def dist(x,y):

    return sum( (x[i]-y[i])*(x[i]-y[i]) for i in xrange(len(y)) )**.5
    

def getWeightedList2(pt,errs):
    '''
    returns a dictionary of distances->values
            and a sorted dictionary of distances
            so d[distDict[1]] is the value of the nearest valid point
    '''
    distList = []
    for ept in errs:
        d = dist(pt,[ept[0],ept[1]])
        distList.append([d,ept[2]])
    return sorted(distList)
    
def matInterpolator2(ErrHash, ranges,res):
    '''
    Interpolates missing data in a heat map by calculating the weighted average
    of the 4 nearest available points, using newer and statistically more accepted
    shephard's method
    '''

    x,y = [int(ranges*res),int(ranges*res)]
    print x,y
    errs = []
    mat = []
    for i in xrange(x):
        mat.append([0.0]*y)
        
    for conc in ErrHash:
        errs.append([conc[0],conc[1],ErrHash[conc]])

    import sys
    
    print '|'+' '*40 + ' |'
    print '|',
    printstep = (x*y)/20.0

    print errs
    #weighted average the top 4 points
    for i in xrange(y):
        for j in xrange(x):
            
            if i*j%int(printstep) == 0: 
                print '=',
                sys.stdout.flush()

  
            #sD is a sorted list [[distance0,error0]...[distanceN,errorN]]
            sD = getWeightedList2([i/res,j/res],errs)
            
            #shepard's method 1968, modified to only consider nearest numpts
            numpts = len(sD)
            p = 1.5
            weightedSum = 0.0
            #error == 0, means this point is definitely known, also we cant compute 0^-n as it would be 1/0
            if sD[0][0] == 0.0:weightedSum = sD[0][1]
            else:
                wibottom = 0.0
                for l in range(numpts):
                        wibottom = wibottom+sD[l][0]**-p
                for k in range(0,numpts): # current just use the nearest 4
                    weightedSum = weightedSum + sD[k][0]**-p/wibottom * sD[k][1]
            mat[i][j] = weightedSum
    return mat



def heatMap(avgErrorMap, ranges,titlestr="",maxError = -1.0,fig=0):
    
    axis = array([i/res for i in range(len(avgErrorMap))])

    
    pcolormesh(axis,axis,avgErrorMap)
    ylabel("X")
    xlabel("Y")
    
    colorbar()
    title(titlestr)
    

    #savefig("/home/lee/Desktop/"+title+".png")
    #clf()

def randomClusters():
    part= 100
    clu = 10
    d = 2
    pts=[]
    for i in range(clu):
        variance = randn()*randn()*.1
        means = [(random()-.5)*6 for b in xrange(d)] 
        for j in range(part):
            p =[0.0]*d
            for k in xrange(d):
                p[k]=(randn()*variance+means[k])
            pts.append(p)
    return pts


def decode(p):
    phat = (p[0]*4+p[1]*4)
    return int(p[0]+.5),int(p[1]+.5)
    
def projectAndDecode(dat,R):
    decodedpts = []
    minimum = 100.0
    maximum = -100.0
    for r in dat:

        dec = decode(r)
        if dec[0]<minimum:minimum = dec[0]
        if dec[1]<minimum:minimum = dec[1]
        if dec[0]>maximum:maximum = dec[0]
        if dec[1]>maximum:maximum = dec[1]
        decodedpts.append(dec)
        
        
    return decodedpts,minimum,maximum
  

def createBuckets(projPoints,mn):
    bucketCounts = {}
    for p in projPoints:
         p = (p[0]+ (0-mn),p[1]+(0-mn))
         if bucketCounts.has_key(p):
            bucketCounts[p] = bucketCounts[p]+1
         else:
            bucketCounts[p] = 1 
    return bucketCounts  
    





print "Generating Data"
points = randomClusters()


R=[[1.,0.],[0.,1.]]
projPoints,mn,mx = projectAndDecode(points,R)
bucketCounts = createBuckets(projPoints,mn)
res = 50.

rang = mx+(0-mn) 

print "interpolating data"
bucketCountsInterp = matInterpolator2(bucketCounts,float(rang),res)
figure(0)
heatMap(array(bucketCountsInterp), res)
figure(1)
scatter([p[1]-mn for p in points],[p[0]-mn for p in points])
show()


