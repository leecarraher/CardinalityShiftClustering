#gcc -c -fPIC leechDecoder.c
#gcc -shared leechDecoder.o -o leechDecoder.so


from ctypes import *
import os

leechdec = cdll.LoadLibrary(os.getcwd()+'/leechDecoder.so')


leechdec.decode.argtypes= [POINTER(c_float),POINTER(c_float)]
leechdec.decode.restype = c_ulong

#leechdec.decode2.argtypes= [POINTER(c_float)]
#leechdec.decode2.restype = c_int

floatArray =c_float*24
def algo2(pts): 

    cx=floatArray(*pts)
    d = c_float()
    result = leechdec.decode(cx,byref(d))

    return long(result)

#def algo1(pts):
#    cx=floatArray(*pts)
#    result = leechdec.decode2(cx)
#    return int(result)


def flatten(l):
    out = []
    for s in l:
        out.extend(s)
    return out
    
    
    
    
    
    
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
    
    #our range is 0-2, so 2 bits each
    if hrd>nrd:
        s = 0
        for i in range(8):
            if nr[i]>2 or nr[i]<0:nr[i]=0
            s = s+((1<<(i*2)) *((nr[i]*2)%4))
    else:
        s = 0
        for i in range(8):
            if nr[i]>2 or hr[i]<0:hr[i]=0
            s = s+((1<<(i*2)) *((hr[i]*2)%4))
    return s



    
    
    
    
    
    
    


def smoothList(list,strippedXs=False,degree=25):  
     if strippedXs==True: return Xs[0:-(len(list)-(len(list)-degree+1))]  
     smoothed=[0]*(len(list)-degree+1)  
     for i in range(len(smoothed)):  
         smoothed[i]=sum(list[i:i+degree])/float(degree)
     aug = [list[i] for i in range(degree-1)]
     avg = sum(aug)/float(len(aug))
     aug = [aug[i]/avg for i in range(len(aug))]
     return aug+ smoothed





def testLatticeDecoders(l,dec):
    import random
    import leechencoder
    import leechdecoder
    errors = []
    weights = [0 for i in range(25)]
    decErrors = [0,0,0,0]
    distances = []
    dat = []
    for e in range(1,800):
        c = 0
        avgdist = 0.0
        for t in range(0,l):
            toterr = 0
            dataBitVec = random.randint(0,2**24-1)
            data = [int(i) for i in bin(dataBitVec).split('b')[1]]
            while len(data)<24:data = [0]+data #pad it out
            #this is a lattice point center
            sent= leechencoder.leechEncoder(data)       
            received =[]
            
            pointErrors = random.randint(1,23)
            overallDistance = (float(e)/30.0)/float(pointErrors)
            errlst = []

            for i in range(12):
                received.append(sent[i])            
            cleanCW = dec(flatten(received))

            # some random errors
            for i in range(pointErrors):
                pl = random.randint(0,23)
                while errlst.__contains__(pl):pl = random.randint(0,23)
                errlst.append(pl)
                if pl%2:
                    received[int(pl/2.0)] = (received[int(pl/2.0)][0]+overallDistance ,received[int(pl/2.0)][1] )
                else:
                    received[int(pl/2.0)] = (received[int(pl/2.0)][0] ,received[int(pl/2.0)][1]-overallDistance )

            codeword = dec(flatten(received))
            weight = 0.0
            winningDecoder = 1
            
            cwstring = [int(i) for i in bin(codeword).split('b')[1]]
            
            weights[sum(cwstring)]=weights[sum(cwstring)]+1

            dist = sum([(received[i][0]-sent[i][0])**2+(received[i][1]-sent[i][1])**2 
            for i in range(12)])
            avgdist = avgdist + dist/float(l)
            if cleanCW-codeword:
                c = c+1

        if e%80==0:print avgdist,c
        dat.append((avgdist,1.0-c/float(l)))
        distances.append(avgdist)
        errors.append(1.0-c/float(l))
    return dat



        #0000 1000 0100 1100   0010 1010 0110 1110
        #0001 1001 0101 1101   0011 1011 0111 1111
pts16 = [(1.0,7.0),(3.0,1.0),(3.0,5.0),(5.0,7.0),  (5.0,3.0),(7.0,5.0),(7.0,1.0),(1.0,3.0),
       (3.0,7.0),(5.0,1.0),(5.0,5.0),(7.0,7.0),(7.0,3.0),(1.0,5.0),(1.0,1.0),(3.0,3.0)]

pts16 = [(pt[0]/8.0,pt[1]/8.0) for pt in pts16]#double it
#pts16 = [(pts16[i][0]/5.6569,pts16[i][1]/5.6569) for i in range(16)]     


def distance(cp,pt):
    '''
    simple euclidean distance formula, can be changed if noise is not gaussian
    '''
    return ((cp[0]-pt[0])**2.0 + (cp[1]-pt[1])**2.0)**.5


def QAM16Dec(r):
    ret = []
    for i in range(6):
        mindist = 1000.0
        minpt = 0
        it = 0
        for pt in pts16:
            
            d = distance((r[i*2],r[i*2+1]),pt)
            if d<mindist:
                mindist = d
                minpt = it
            it = it +1
        cw = [int(c) for c in bin(minpt).split('b')[1]]
        cw = [0]*(4-len(cw))+cw
        cw.reverse()
        ret.extend( cw  )
    return ret


def QAM16Enc(r):
    ret = []
    for i in range(6):
        #last one is *8
        #next is *4
        s = r[i*4+3]*8+r[i*4+2]*4+r[i*4+1]*2+r[i*4]
        ret.append(pts16[s])
    return ret



pts64 = [(.5,3.5),(1.5,.5),(1.5,2.5),(2.5,3.5),(2.5,1.5),(3.5,2.5),
(3.5,.5),(.5,1.5),(1.5,3.5),(2.5,.5),(2.5,2.5),(3.5,3.5),(3.5,1.5),
(.5,2.5),(.5,.5),(1.5,1.5),#00
        (-3.5,3.5),(-2.5,.5),(-2.5,2.5),(-1.5,3.5),(-1.5,1.5),(-.5,2.5),
        (-.5,.5) ,(-3.5,1.5),(-2.5,3.5),(-1.5,.5),(-1.5,2.5),(-.5,3.5),
        (-3.5,2.5),(-3.5,.5),(-.5,1.5),(-2.5,1.5),#10
        (-3.5,-.5),(-2.5,-3.5),(-2.5,-1.5),(-1.5,-.5),(-1.5,-2.5),
        (-.5,-1.5),(-.5,-3.5),(-3.5,-2.5),(-2.5,-.5),(-1.5,-3.5),(-1.5,-1.5),
        (-.5,-.5),(-.5,-2.5),(-3.5,-1.5),(-3.5,-3.5),(-2.5,-2.5),#01
        (.5,-.5),(1.5,-3.5),(1.5,-1.5),(2.5,-.5),(2.5,-2.5),(3.5,-1.5),
        (3.5,-3.5),(.5,-2.5),(1.5,-.5),(2.5,-3.5),(2.5,-1.5),(3.5,-.5),
        (3.5,-2.5),(.5,-1.5),(.5,-3.5),(1.5,-2.5) ]#11

pts64 = [(pt[0]/8.0,pt[1]/8.0) for pt in pts64]#double it


def QAM64Dec(r):
    ret = []
    for i in range(4):
        mindist = 1000.0
        minpt = 0
        it = 0
        for pt in pts64:
            d = distance((r[i*2],r[i*2+1]),pt)
            if d<mindist:
                mindist = d
                minpt = it
            it = it +1
        cw = [int(c) for c in bin(minpt).split('b')[1]]
        cw = [0]*(6-len(cw)) + cw# zero pad
        cw.reverse()
        ret.extend( cw  )
    return ret 


def ebN0Leech(data,dec,EbN0dB ):
    from random import normalvariate,randint
    from math import log10,log
    from leechencoder import leechEncoder
    from pylab import arange
    

    M=16.0 #Number of Constellation points M=2^k
    Rm=log(M,2.0) # 
    Rc=1.0/2.0  #Rc = code rate for a coded system. may also be

    BER = [0.0]*len(EbN0dB)# BER values for each Eb/N0

    index = 0
    for i in EbN0dB:
        #Channel Noise for various Eb/N0
        #Adding noise with variance according to the required Eb/N0
        EbN0 = 10.0**(float(i)/10.0); #Converting Eb/N0 dB value to linear scale
        noiseSigma =(1./(2*Rm*Rc*EbN0))**.5 #Standard deviation for AWGN Noise
        
        #Creating noise for adding with coded signal
        totalBER = 0
        for j in range(0,N/24):
            sent= leechEncoder(data[j*24:(j+1)*24])
            sent = flatten(sent)
            y = dec(sent)
            received = [sent[k]+8*noiseSigma*normalvariate(0.0,1.0) for k in range(24)]
            x = dec(received)
            
            totalBER = totalBER + sum([int(k) for k in bin(x^y).split('b')[1]])

        res = float(totalBER)/ (float(len(data)*2))
        if res !=0.0:
            BER[index] = log10(res)
        else:
            BER[index] = min(BER)

        print BER[index],i
        index = index+1
    return BER
    
    
def ebN0E8(data,EbN0dB):
    from random import normalvariate,randint
    from math import log10,log
    from pylab import arange
    

    M=2.0 #Number of Constellation points M=2^k
    Rm=log(M,2.0) # 
    Rc=1.0/2.0  #Rc = code rate for a coded system. may also be

    BER = [0.0]*len(EbN0dB)# BER values for each Eb/N0

    index = 0
    for i in EbN0dB:
        #Channel Noise for various Eb/N0
        #Adding noise with variance according to the required Eb/N0
        EbN0 = 10.0**(float(i)/10.0); #Converting Eb/N0 dB value to linear scale
        noiseSigma =(1./(2*Rm*Rc*EbN0))**.5 #Standard deviation for AWGN Noise
        
        #Creating noise for adding with coded signal
        totalBER = 0
        for j in range(0,N/8):
            scaledDat = [0,0,0,0,0,0,0,0]
            received = [0,0,0,0,0,0,0,0]
            y = 0
            for k in range(8):
                scaledDat[k] = 2.0*(data[j*8+k])
                y = y+((1<<(k*2)) *((scaledDat[k]*2)%4))
                received[k] = scaledDat[k]+noiseSigma*normalvariate(0.0,1.0)
            #print received,scaledDat
                
                
            x = int(decodeE8(received))
            y = int(y)
            totalBER = totalBER + sum([int(k) for k in bin(x^y).split('b')[1]])

        res = float(totalBER)/ float(len(data))
        if res !=0.0:
            BER[index] = log10(res)
        else:
            BER[index] = min(BER)

        print BER[index],i
        index = index+1
    return BER
    
    
    
    
    
    
    
    

def ebN0QAM(N,M,QAMEnc,QAMDec,EbN0dB ):
    from random import normalvariate,randint
    from math import log10,log
    from pylab import arange

    data=  [randint(0,1) for i in range(N) ] #random data

    Rm=log(M,2.0)
    Rc=1.0  #Rc = code rate

    BER = [0.0]*len(EbN0dB)# BER values for each Eb/N0

    index = 0
    for i in EbN0dB:
        #Channel Noise for various Eb/N0
        #Adding noise with variance according to the required Eb/N0
        EbN0 = 10.0**(float(i)/10.0); #Converting Eb/N0 dB value to linear scale
        noiseSigma =(1./(2*Rm*Rc*EbN0))**.5 #(1./(2.0*Rm*Rc*EbN0))**.5 #Standard deviation for AWGN Noise
        
        #Creating noise for adding with coded signal
        totalBER = 0
        for j in range(0,N/24):
            datatemp = data[j*24:(j+1)*24]
            sent = flatten(QAMEnc(datatemp))
            received = [sent[k]+noiseSigma*normalvariate(0.0,1.0) for k in xrange(len(sent))]
            
            decoded = QAMDec(received)   
            totalBER = totalBER + sum([decoded[k]^datatemp[k] for k in xrange(24)])

        res = float(totalBER)/ float(len(data))
        if res !=0.0:
            BER[index] = log10(res)
        else:
            BER[index] = min(BER)

        print BER[index],i
        index = index+1
    return BER
def theoreticalQAM(M,EbN0dB):
    '''
    this function generates the theoretical error rate for a QAM modulated
    signal where the number of lattice points is M and of the form
    2^k=M
    '''
    from scipy import special
    k = 1/((2/3.0)*(M-1))**.5
    return [log(2*(1-1/M**(.5))*special.erfc(k*(10**(eb/10))**.5) - 
    (1-2/M**.5 + 1.0/M)*(special.erfc(k*(10**(eb/10))**.5))**2)/log(10) for eb in EbN0dB]
    
def ebN0Golay(N,EbN0dB):
    from random import normalvariate,randint
    from math import log10,log
    from pylab import arange
    from golay import decodepts,encodepts

    
    data=  [randint(0,1) for i in range(N) ] #random data

    M=16.0 #Number of Constellation points M=2^k
    Rm=log(M,2.0) # 
    Rc=1.0/2.0  #Rc code rate

    BER = [0.0]*len(EbN0dB)# BER values for each Eb/N0

    index = 0
    for i in EbN0dB:
        #Channel Noise for various Eb/N0
        #Adding noise with variance according to the required Eb/N0
        EbN0 = 10.0**(float(i)/10.0); #Converting Eb/N0 dB value to linear scale
        noiseSigma =(1./(2*Rm*Rc*EbN0))**.5 #Standard deviation for AWGN Noise
        
        #Creating noise for adding with coded signal
        totalBER = 0
        for j in range(0,N/12):
            datatemp = data[j*12:(j+1)*12]
            sent = encodepts(datatemp)
            received = [(sent[k][0]+noiseSigma*normalvariate(0.0,1.0),
            sent[k][1]+noiseSigma*normalvariate(0.0,1.0) ) for k in xrange(6)]           
            decoded = decodepts(received)
            totalBER = totalBER + sum([decoded[k]^datatemp[k] for k in xrange(12)])

        res = float(totalBER)/ float(len(data))
        if res !=0.0:
            BER[index] = log10(res)
        else:
            BER[index] = min(BER)
        
        
        print BER[index],i
        index = index+1
    return BER

if __name__ == '__main__':
    pass
    
    from random import random, randint
    import leechdecoder
    import leechencoder
    import pylab
    
    '''dat = testLatticeDecoders(1000,algo1)
    dat.sort()
    pylab.plot([ d[0] for d in dat],[ d[1] for d in dat],label="Smoothed Bounded Distance Algo. 1")


    dat = testLatticeDecoders(1000,algo2)
    dat.sort()
    pylab.plot([ d[0] for d in dat],[ d[1] for d in dat],label="Smoothed Bounded Distance Algo. 2")

    pylab.legend(loc='upper right')
    pylab.xlabel('Average Distance')
    pylab.ylabel('Collisions Probability')

    pylab.show()
    
    
    '''
    N= 384000/4
    data=  [randint(0,1) for i in range(N) ] #random data

    EbN0dB = pylab.arange(-1.0,17.3,.7)
    
    
    #BERQ64 = ebN0QAM(N,64,QAM64Enc,QAM64Dec,EbN0dB)
    #BERQ16 = ebN0QAM(N,16,QAM16Enc,QAM16Dec,EbN0dB)
    
    #BERLE8 = ebN0E8(data,EbN0dB)
    BERL =ebN0Leech(data,algo2,EbN0dB)
    #BERL2 = ebN0Leech(data,algo2,EbN0dB)
    #BERGol = ebN0Golay(N,EbN0dB)
    


    pylab.plot(EbN0dB, BERL,'v-',label='Leech hex -.2db and QAM16')
    #pylab.plot(EbN0dB, BERL2,label='Leech hex algo2 and QAM16')
    #pylab.plot(EbN0dB, BERLE8,'x-',label='E8 2PSK')
    #pylab.plot(EbN0dB, BERGol,'v',label='Golay Code and QAM16')
    #pylab.plot(EbN0dB, BERQ64,'o-',label='Unencoded QAM64')
   #pylab.plot(EbN0dB, BERQ16,'D-',label='Unencoded QAM16')
    pylab.grid(True)
    pylab.title("Error Correcting Performance ("+str(N)+" Samples)")
    pylab.legend(loc='lower left')
    pylab.xlabel('EB/NO (dB)')
    pylab.ylabel('Probability of Bit Error - log10(Pb)')
    pylab.show()
    


