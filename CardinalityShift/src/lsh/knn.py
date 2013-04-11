'''
Copyright 2010 Lee Carraher. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Lee Carraher ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Lee Carraher OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of Lee Carraher.


author: lee carraher

this is an ad-hoc version of the e2lsh (by Alexandr Andoni and 
Piotr Indyk) program in python it uses the leech decoder compiled
obj for the lsh function it then concatenates the hashes of the 
vectors for a given object at 240 bits and stores this with a 
pointer bach to the originating object in a hash table.

calling this module as main will build the hash table for items 
in the pgm directory and store it as storedDatabase.pickle

importing this module allows a stored database to be used for 
imaging matching unseen or test data, without having to 
regenerate the database each time, ie for 1 off tests or 
something like that. 

'''

#import shelve
from pylab import dot,array,log,histogram,plot,show,eye

# leech library stuff
from ctypes import cdll,c_ulong,c_float,POINTER,byref

import time,random,numpy,os

leechdec = cdll.LoadLibrary(os.getcwd()+'/leechDecoder.so')
#siftmatch = cdll.LoadLibrary(os.getcwd()+'/match.so')
#siftmatch.FindMatches2.argtypes= [POINTER(c_char),POINTER(c_char)]
#siftmatch.FindMatches2.restype = c_int



leechdec.decode.argtypes= [POINTER(c_float),POINTER(c_float)]
leechdec.decode.restype = c_ulong


def distance(x,y,d):
    dist =0.0
    for i in xrange(d):
        dist = dist+(x[i]-y[i])*(x[i]-y[i])
    return dist**.5
    

def linearSearchNN(q,DB,d,minVal=0.0,maxVal=1.0):
    minDist = 1000000.0
    
    for i in xrange(len(DB)):
        dist = distance(q,DB[i],d)
        if dist<minDist:
            minDist = dist
            minq = i
    return minq



floatArray =c_float*24
l = 2
tblSize = 4096**l



def cnt(v):
    c =  ((v & 0xfff) * 0x1001001001001 & 0x84210842108421) % 0x1f
    c += (((v & 0xfff000) >> 12) * 0x1001001001001 & 0x84210842108421)%0x1f
    return c
    
maxlen =0
'''
maybe a cleverway to partition hashes
using what we know about octads
2576^n options
2*759^n
2*2^n
and all permutations
'''    
'''def ELFHash(key):
    #key=str(key)
    hash = 0
    x = 0
    while key>0:
        part = key&((1<<24)-1)
        key = key>>24
        
        t = cnt(part)
        #largest subset is 12
        if t==12:
            hash = (hash<<12)+(part%4079)
        elif t ==0 or t ==24:
            hash = (hash<<10)+(part%751)
        else:
            hash = (hash<<1)+(part%2)

    return hash
'''

def ELFHash(key):
    #dist = key[1] 
    #key = key[0]
    key = str(key)
    h = 0
    x    = 0
    for i in range(len(key)):
        h = (h << 4) + ord(key[i])
        x = h & 0xF0000000
        if x != 0:
            h ^= (x >> 24)
        h &= ~x
    return h%tblSize#,dist


def scale(d,low,high,newLow,newHigh):
    k = (newHigh - newLow) / (high - low)
    for i in xrange(len(d)):        
        d[i] = newLow+(d[i]-low)*k



def decode(pts):

    scale(pts,0.,1.,0.0,8.0)
    cx=floatArray(*pts)
    d = c_float()
    result = leechdec.decode(cx,byref(d))
    r = long(result) & 0xffffff
    if not (cnt(r)==8 or cnt(r)==12 or cnt(r)==16 or cnt(r)==24) :  print long(result)
    return long(result),d.value
    
#def match(file1,file2):
#    f1=(c_char*len(file1))(*file1) 
#    f2=(c_char*len(file2))(*file2) 
#    result = siftmatch.FindMatches2(f1,f2)
#    return int(result)
    
    
    

def dist(r,d):
    tot = 0.0
    for i in range(128):
        tot = tot+(r[i]-d[i])**2
    return tot**.5
    
    

def hashcompare(v1,v2):
    incommon = 0
    for i in v1:

        for j in v2:
            if i ==j:incommon=incommon+1

    return incommon
    


def buildDB(data,minVal=0.0,maxVal=1.0):
    n = len(data)
    d = len(data[0])
    #parameters
    P1=.01779
    P2=.0000156
    rho = log(P1)/log(P2)
    k =  int((log(n)/log( 1/P2)+.5))
    sim = int(log(n)+.5)
    
    M = generateRandomProj(k,sim,d)
    
    print n**rho,k,sim
    
    DBs = [{} for s in xrange(sim)]
    for s in xrange(sim):
        for i in xrange(n):
            hashVal = decodeGt24(data[i],M[s],minVal,maxVal)
            if DBs[s].has_key(hashVal):
                DBs[s][hashVal].append(i)
            else:
                DBs[s][hashVal]=[i]
                    
    return DBs,M
    
def query(q,DBs,Ms,n,l=1000,minVal=0.0,maxVal=1.0):
    #parameters
    P1=.01779
    P2=.0000156
    rho = log(P1)/log(P2)
    sims = len(Ms)#number of multi-runs for whp
    k = len(Ms[0])#number of random matrix projections per run
    
    candidates = set()
    #first iterate over the runs
    for s in xrange(sims):
        #next iterate over the n^rho nearby points
        hashVal = decodeGt24(q,Ms[s],minVal,maxVal)
        if DBs[s].has_key(hashVal):
        
            for c in DBs[s][hashVal]:
                candidates.add(c)
                
        for r in xrange(int(n**rho+.5)):
            hashVal = decodeGt24(q,Ms[s],minVal,maxVal,True)
            
            if DBs[s].has_key(hashVal):
                for c in DBs[s][hashVal]:
                    candidates.add(c)
            if len(candidates)>2*l:return candidates
            
    return candidates


    
def flatten(l):

    out = []
    for s in l:
        out.extend(s)
    return out

#end leech library stuff
#constant scaling factor 1/sqrt(24)
sc = 1.0/((24)**.5)
maxint = 0

def convert2Bytes(b):
    '''
    shelves and dbs only hold strings and ints, we could naively 
    create strings of the longs, but then it assumes each is 8 bits
    this conversion is more efficient for storage space
    '''
    b = bin(b)[2:]
    ret = ''
    for i in range(len(b)/8):
        ret = ret + chr(int(b[i*8:(i+1)*8],2))
    return ret
    
    

def generateRandomProj(k,sim,d):
    M = [[] for i in xrange(sim)]
    for i in range(sim):
        for j in range(k):
            M[i].append( array([[gauss(0,1)*sc for a in xrange(24)] for b in xrange(d)]))
    return M
    
    
def decodeD8(r):
    rInt = [0,0,0,0,0,0,0,0]
    rDist = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    for i in range(8):
        tr = int(r[i])
        if r[i]>tr+.5:
            rInt[i]=tr+1
            rDist[i] = (tr+1)-r[i]
        else:
            rInt[i]=tr
            rDist[i] = r[i]-tr
    if sum(rInt)%2==0:
        return rInt
    else:
        rDist.index(max(rDist))
        if rInt[i]>r[i]:
            rInt[i]=rInt[i]-1
        else: 
            rInt[i]=rInt[i]+1
        
    return rInt
    
def decodeE8(r):
    '''
        this decoder uses the D8 partition to decode E8 lattice
        it returns an integer label for the lattice point 
    '''


    nr = decodeD8(r)
    hr = [k+.5 for k in decodeD8([j-.5 for j in r])]
    
    nrd = sum([(nr[i]-r[i])**2 for i in range(8)])
    hrd = sum([(hr[i]-r[i])**2 for i in range(8)])
    
    #our range is 0-2, so normalized to 2
    if hrd>nrd:
        s = 0
        for i in range(8):
            s = s+((1<<(i*2)) *((nr[i]*2)%4))
        return s
    else:
        s = 0
        for i in range(8):
            s = s+((1<<(i*2)) *((hr[i]*2)%4))
        return s
        
    

def decodeGt24(v,M,minVal,maxVal,rand=False):
    scale(v,minVal,maxVal,0.0,1.0)
    #distance = 0.0
    hashVal = 0
    i=0
    while i < len(M):
        vproj = dot(v,M[i]).tolist()
        #if rand:vproj=[vhat+gauss(0,.047)/sc for vhat in vproj]
        dec = decode(vproj)
        #distance = distance + 1.0/(dec[1]+1.0)
        hashVal=(hashVal<<37) + dec[0]
        i = i+1
        
    return ELFHash(hashVal)#onvert2Bytes(hashVal)


#<-----------Testing stuff--------------->    


def decodeE8gt8(v,M,hi,low):
    scale(v,hi,low,0.0,2.0)
    hashVal = 0
    i=0
    while i < len(M):
        vproj = dot(v,M[i]).tolist()
        dec = int(decodeE8(vproj))
        #distance = distance + 1.0/(dec[1]+1.0)
        hashVal=hashVal + (dec<< (16*i) )
        i = i+1
        
    return ELFHash(hashVal)
 

from random import gauss,random

def testCollisionsE8(n,d=24):
    
    M =[array([[gauss(0,1)/(24**.5) for a in range(24)] for b in range(d)]) for k in range(12)]
    S = [0.0]*n
    C = [0]*n
    #generate distances and buckets
    i=0
    while i <n:
        p = [random() for j in xrange(d)]
        q = [p[j] + (gauss(0,1)/(d**.5)) for j in xrange(d)]
        
        S[i]=distance(p,q,d)
        j=0
        m = 3

        while j < len(M)/m:
            C[i]= C[i] or decodeE8gt8(p,M[j*m:(j+1)*m],0.0,1.0)==decodeE8gt8(q,M[j*m:(j+1)*m],0.0,1.0)
            j = j +1
        i=i+1
    
    ranges = histogram(S,25)[1]   
    bucketsCol = [0]*len(ranges)
    bucketsDis = [0]*len(ranges)


    #fill buckets with counts 
    for i in xrange(n):
        k = len(ranges)-1
        while S[i] < ranges[k]:k=k-1
        if C[i]:bucketsCol[k]=bucketsCol[k]+1
        bucketsDis[k] = bucketsDis[k]+1
        
    print bucketsDis
    print bucketsCol
    print ranges

    p = [0]*len(ranges)
    for m in range(len(ranges)):
        if bucketsDis[m]>0:
            p[m]=bucketsCol[m]/float(bucketsDis[m])

    plot(ranges,p) 

def gaussInv(mu,sig):
    from random import gauss,random
    s = 1.
    if random()>.5: s = -1. 
    if random()>.5:
        #lhs
        g = gauss(-mu/2.,sig)
        while g < -mu:g = gauss(-mu/2.,sig)
        return g*s
    #rhs
    g = gauss(mu,sig)
    while g > mu/2.:g = gauss(mu/2.,sig)
    return g*s


def testCollisionPr(n ,d=24):
    M =[array([[gauss(0,1)/(24**.5) for a in range(24)] for b in range(d)]) for k in range(4)]
    #M =[eye(24)]#[array([[gauss(0,1)/(24**.5) for a in range(24)] for b in range(d)])]
    S = [0.0]*n
    C = [0]*n
    print "data generated"
    #generate distances and buckets
    i=0

    while i <n:
       
        p = array([(random()) for j in xrange(d)])
        
        q =array( [p[j] + (gaussInv(3.0,1)/(d**.5)) for j in xrange(d)])
        #q =array( [p[j] + (gauss(0.0,1.0)/(d**.5)) for j in xrange(d)])
        #q =array( [p[j] + (( random() -.5)/(d**.5)) for j in xrange(d)])
        S[i]=distance(p,q,d)
        j=0
        m = 2

        while j < len(M)/m:
            C[i] = C[i] or decodeGt24(p,M[j*m:(j+1)*m],0.0,1.0)==decodeGt24(q,M[j*m:(j+1)*m],-0.0,1.0)
            j = j +1
            
        i = i+1
        
    ranges = histogram(S,25)[1]   
    bucketsCol = [0]*len(ranges)
    bucketsDis = [0]*len(ranges)


    #fill buckets with counts 
    
    for i in xrange(n):
        k = len(ranges)-1
        while S[i] < ranges[k]:k=k-1
        if C[i]:bucketsCol[k]=bucketsCol[k]+1
        bucketsDis[k] = bucketsDis[k]+1
        
    print bucketsDis+bucketsCol
    print ranges

    p = [0]*len(ranges)
    for m in range(len(ranges)):
        if bucketsDis[m]>0:
            p[m]=bucketsCol[m]/float(bucketsDis[m])

    plot(ranges,p) 


if __name__ == '__main__':
    pass
    #testCollisionsE8(100000)
    testCollisionPr (100000)
    show()
    

