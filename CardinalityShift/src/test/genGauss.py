
from numpy import random

def getDataPoints(part,d,clu):
    ret = [[] for j in xrange(part*clu)]
    clusterCenters = []
    for i in range(clu):
        variance = .5+random.random()*4
        means = [random.random(1)*10 for b in range(d)] 
        clusterCenters.append(means)
        for j in range(part):
            p = [0]*d
            for k in xrange(d):
                pts = random.random()
                p[k]= (pts*variance+means[k])[0]
            ret[i*part+j]=p
    return ret,clusterCenters

if __name__ == '__main__':
    pass

    part= 3
    clu = 5
    d = 8
    ret,cluCntrs = getDataPoints(part,d,clu)
    
    print "x y z label"
    for i in range(clu):
        for j in range(part):
            for k in xrange(d):
                print ret[i*part+j][k],
            print chr(i+97)