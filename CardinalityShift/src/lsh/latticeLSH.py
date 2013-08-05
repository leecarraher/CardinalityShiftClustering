'''
Created on May 7, 2013

@author: lee
'''
from math import log,sqrt,ceil
from numpy import array,dot
import os
from random import gauss,random,randint

from ctypes import *

leechdec = cdll.LoadLibrary(os.getcwd()+'/leechDecoder.so')
leechdec.decode.argtypes= [POINTER(c_float),POINTER(c_float)]
leechdec.decode.restype = c_ulonglong

floatArray =c_float*24


numberofdataPoints = 10000
numberofclusers = 20
dimensions = 100


def getDataPoints(part,d,clu):
    ret = [[] for j in xrange(part*clu)]
    clusterCenters = []
    for i in range(clu):
        variance = .5+random()*2
        means = [random()*8 for b in range(d)] 
        clusterCenters.append(means)
        for j in range(part):
            p = [0]*d
            for k in xrange(d):
                pts = random()
                p[k]= (pts*variance+means[k])
            ret[i*part+j]=p
    return ret,clusterCenters




def decode(pts):
    cx=floatArray(*pts)
    d = c_float()
    result = leechdec.decode(cx,byref(d))
    return int(result),d.value
    
P1=.01779
P2=.0000156
rho = log(P1)/log(P2)
sc = 1.0/(24**.5)


from pylab import array
def decodeGt24(v,k = 1):

    distance = 0.0
    hashVal = 0


    for i in xrange(k):
        #print v
        R = array([[ gauss(0,sc)/7 for a in xrange(24)] for b in xrange(len(v))])
        
        vproj=dot(v,R).tolist()
        vproj = [vv+4. for vv in vproj]
        #print max(vproj),min(vproj)
        
        dec = decode(vproj)
        distance = distance + 1.0/(dec[1]+1.0)
        #xor ignors is directionless for concatonations
        hashVal=hashVal^dec[0]
        
        
    return hashVal



def buildDB(data):
    n = len(data)

    k =  int(ceil((log(n)/log( 1.0/P2)+.5)))# length of concatenations, increases specificity

    DB = {}
    for i in xrange(n):
        hashVal = decodeGt24(data[i],k)
        if DB.has_key(hashVal):
            DB[hashVal].append(i)
        else:
            DB[hashVal]=[i]

    return DB
    

def query(q,DBs,n,l=10,k=1):
        #parameters

    s = len(DBs)
    candidates = set()
    
    for r in xrange(int(n**rho+.5)):
            hashVal = decodeGt24(q,k)
            if DBs[s].has_key(hashVal):
                
                for c in DBs[s][hashVal]:
                    candidates.add(c)
            if len(candidates)>2*l:return candidates  
            
            
    return candidates





if __name__ == '__main__':
    pass

    data,cntrs = getDataPoints(numberofdataPoints,dimensions,numberofclusers)
    print "building DB"
    DB =  buildDB(data)
    print "Querying DB size:"+str(len(DB)) + " pts: " + str(len(data))
    
    for i in xrange(len(cntrs)):
        print query(cntrs[i],DB,len(data))


        