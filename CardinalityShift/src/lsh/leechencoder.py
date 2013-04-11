#access codewords by their information vector indeces H6CodeWords[firstchar][secondchar][thirdchar]
#ex. wW1 - > H6CodeWords[2][3][1] = [0,1,0] => the codeword is wW1010
H6CodeWords = [[[[0,0,0],[1,1,1],[2,2,2],[3,3,3]],[[1,2,3],[0,3,2],[3,0,1],[2,1,0]],[[2,3,1],[3,2,0],[0,1,3],[1,0,2]],[[3,1,2],[2,0,3],[1,3,0],[0,2,1]]],[[[1,3,2],[0,2,3],[3,1,0],[2,0,1]],[[0,1,1],[1,0,0],[2,3,3],[3,2,2]],[[3,0,3],[2,1,2],[1,2,1],[0,3,0]],[[2,2,0],[3,3,1],[0,0,2],[1,1,3]]],[[[2,1,3],[3,0,2],[0,3,1],[1,2,0]],[[3,3,0],[2,2,1],[1,1,2],[0,0,3]],[[0,2,2],[1,3,3],[2,0,0],[3,1,1]],[[1,0,1],[0,1,0],[3,2,3],[2,3,2]]],[[[3,2,1],[2,3,0],[1,0,3],[0,1,2]],[[2,0,2],[3,1,3],[0,2,0],[1,3,1]],[[1,1,0],[0,0,1],[3,3,2],[2,2,3]],[[0,3,3],[1,2,2],[2,1,1],[3,0,0]]]]

H6CodeWordsRev = [[[[0,0,0],[1,1,1],[2,2,2] ,[3,3,3]],[[2,3,1],[3,2,0] ,[0,1,3] ,[1,0,2]],[[3,1,2],[2,0,3] ,[1,3,0] ,[0,2,1]],[[1,2,3],[0,3,2] ,[3,0,1] ,[2,1,0]]],[[[3,2,1],[2,3,0],[1,0,3] ,[0,1,2]],[[1,1,0],[0,0,1],[3,3,2],[2,2,3]],[[0,3,3],[1,2,2],[2,1,1],[3,0,0]],[[2,0,2],[3,1,3],[0,2,0] ,[1,3,1]]],[[[1,3,2],[0,2,3],[3,1,0] ,[2,0,1]],[[3,0,3],[2,1,2] ,[1,2,1] ,[0,3,0]],[[2,2,0],[3,3,1],[0,0,2] ,[1,1,3]],[[0,1,1],[1,0,0],[2,3,3] ,[3,2,2]]],[[[2,1,3],[3,0,2],[0,3,1] ,[1,2,0]],[[0,2,2],[1,3,3],[2,0,0] ,[3,1,1]],[[1,0,1],[0,1,0],[3,2,3] ,[2,3,2]],[[3,3,0],[2,2,1],[1,1,2] ,[0,0,3]]]]

def test():
    for i in range(4):
        for j in range(4):
            for k in range(4):
                s = leechencoder.H6CodeWords[i][j][k]
                r = leechencoder.H6CodeWordsRev[s[2]][s[1]][s[0]]
                if r[0]!=i and r[1]!=j and r[2]!=k:print "error"


#there are 4 projections for every character,2 odd and 2 even	
#             e1        e2          o1      o2	
projExpn=[[[0,0,0,0],[1,1,1,1],[0,1,1,1],[1,0,0,0]]  #0 0
	,[[0,0,1,1],[1,1,0,0],[0,1,0,0],[1,0,1,1]]  #1 1
	,[[0,1,0,1],[1,0,1,0],[0,0,1,0],[1,1,0,1]],#w 2
	[[0,1,1,0],[1,0,0,1],[0,0,0,1],[1,1,1,0]]]  #W 3

def partitionBin2GF4(v):
    l = len(v)/2
    ret = [0 for i in range(l)]
    for i in range(l):
        ret[i]=2*v[i*2]+v[i*2+1]
    return ret
        
#multilevel hexacode based golay encoder (totally unoptimized)
def multiLevelGolayConstruct(v):
    '''
        this is a golay code encoder that uses the multilevel construction based on the hexacode
    '''
    hexbits = partitionBin2GF4(v[0:6])
    paritySelectBit = v[6]
    hexacodeword = hexbits + H6CodeWords[hexbits[0]][hexbits[1]][hexbits[2]]
    codeword = []
    #remaining pref. bits choose which projection is choosen between complement and non-complement
    preferenceBits = v[7:]
    parity = 0
    for i in range(len(preferenceBits)):
        p = preferenceBits[i]
        codeword.extend(projExpn[hexacodeword[i]][paritySelectBit*2+p])
        parity=parity+projExpn[hexacodeword[i]][paritySelectBit*2+p][0]

    parity = parity%2
    codeword.extend(projExpn[hexacodeword[5]][paritySelectBit*2+(parity != paritySelectBit)])

  
    return codeword


def parityEncoder(v,paritySelector):
    '''
        this is a parity Encoder (12,11) that appends a bit depending on the parity selector control bit's parity
    '''
    tot = 0
    for i in v:
        tot = tot+i
    v.append(int(tot%2!=paritySelector))
    return v

def repititionEncoder(v):
    '''
        repitition encoder specifies the A/B selector bit
    '''
    return [v for i in range(12)]


#ungerboeck's partitioning
#                    A                                          B                                        A/B
#         A0                   A1                     B0                   B1                            firstbit
#   A00       A11        A10       A01          B00       B11        B10       B01                       secondbit
#A000 A001 A110 A111  A100 A101 A010 A011    B000 B001 B110 B111  B100 B101 B010 B011                    k-parity

#not rotated
#0   1   2    3                                     4    5    6  7
#000 010 100 110                                   001  011 101 111
APts = [(1.0,7.0),(3.0,5.0),(3.0,1.0),(5.0,7.0), (5.0,3.0),(7.0,1.0),(7.0,5.0),(1.0,3.0)]
#APts = [(pt[0]/8.0,pt[1]/8.0) for pt in APts]
#000  010 100 110                                 001  011 101 111
BPts = [(3.0,7.0),(5.0,5.0),(5.0,1.0),(7.0,7.0),(7.0,3.0),(1.0,1.0),(1.0,5.0),(3.0,3.0)]
#BPts = [(pt[0]/8.0,pt[1]/8.0) for pt in BPts]



def ungerbockPtSelect(repEnc,golayEnc,paritySel):
    pts = []
    #since its a reqr that all points be either A/B this is pretty redundant
    ptLat = APts
    if repEnc[0]==1:
        ptLat = BPts
    for i in range(12):    
        index = 2*golayEnc[2*i]+golayEnc[2*i+1]+4*paritySel[i] 
        pts.append(ptLat[index])
       
    return pts

def leechEncoder(v):
    rep = repititionEncoder(v[23])
    gol = multiLevelGolayConstruct(v[0:12])
    par = parityEncoder(v[12:23],v[23])

    ret = ungerbockPtSelect(rep,gol,par)   
    
    return ret

def test():
    import leechDecoder
    import random
    #v=[1,0,1,0,0,0,0,1,1,1,0,1,1,0,0,0,1,1,0,1,1,1,1,0]#Apt_even ,correct position distance and decoding
    #v=[0,1,0,0,1,0,1,1,0,0,1,0,0,0,1,0,1,1,1,0,1,1,1,0]#Apt_odd  ,right position wrong decoding
    #v=[0,0,1,1,0,1,0,1,1,0,1,0,1,1,0,0,0,1,0,1,0,0,1,1]#Bpt_even ,wrong position wrong decoding
    #v=[0,1,0,1,0,0,1,1,0,0,0,0,0,1,1,0,1,0,0,1,1,1,1,1]#Bpt_odd  ,correct position wrong distance,correct decoding
    cc = 0
    for i in range(0,2<<11):
        strWord = bin(i).split('b')[1]
        l = len(strWord)
        s = [0 for i in range(12-l)]
        for j in range(l):
            s.append(int(strWord[j]))
        v=multiLevelGolayConstruct(s)
        
        idx = 0
        if v[23]== 1:idx=2

        tt = multiLevelGolayConstruct(v[0:12])
        if (tt[0]+tt[4]+tt[8]+tt[12]+tt[16]+tt[20]) %2 == 1:idx = idx+1

        if idx>0:
       
            v = [random.randint(0,1) for i in range(24)] 
            trans = leechEncoder(v)
            #deviations by more than 2.46668044 break the correction
            # is this the correct covering radius
            
            trans[1]=(trans[1][0]-.5,trans[1][1]+.1)
            trans[2]=(trans[2][0]-.3,trans[2][1]+.2)
            trans[5]=(trans[5][0]-.9,trans[5][1]+.066689)
            trans[7]=(trans[7][0]-.5,trans[7][1]+.2)
            trans[11]=(trans[11][0]-.1,trans[11][1]+.3)
            trans[3]=(trans[3][0]-.1,trans[3][1]+.2)

            a,b,c,d= leechDecoder.decode(trans)
            #print a
            if b != 0.0:
                
                cw1 = multiLevelGolayConstruct(v[0:12])
                cw2 = c

                if sum([cw1[r]-cw2[r] for r in range(24)]):
                    cc = cc+1
                    print b

    print cc

#64qam sends 2 other binary symbols that select the 1,2,3,4 quadrant
def test2():
    counts = [0,0,0,0]
    cou2s = [0,0,0,0]
    v=[1,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0,1,1,0,1,1,1,1,0]

    for i in range(0,2<<11):
        strWord = bin(i).split('b')[1]
        l = len(strWord)
        s = [0 for i in range(12-l)]
        for j in range(l):
            s.append(int(strWord[j]))
        v=multiLevelGolayConstruct(s)

        idx = 0
        if v[23]== 1:idx=2

        tt = multiLevelGolayConstruct(v[0:12])
        if (tt[0]+tt[4]+tt[8]+tt[12]+tt[16]+tt[20]) %2 == 1:idx = idx+1

        cou2s[idx] = cou2s[idx]+1
        #[1,0,0,1,1,1,1,1,0,0,1,1,0,1,0,1,1,0,0,1,1,1,1,1]
        #W01wW0

        trans = leechEncoder(v)
        import leechDecoder
        a,b,c,d= leechDecoder.decode(trans)

        counts[d]=counts[d]+1
    print counts
    print cou2s
def test3():
    from leechdecoder import decode
    v=[[1,0,1,0,0,0,0,1,1,1,0,1,1,0,0,0,1,1,0,1,1,1,1,0],#Apt_even ,correct position distance and decoding
    [0,1,0,0,1,0,1,1,0,0,1,0,0,0,1,0,1,1,1,0,1,1,1,0],#Apt_odd  ,right position wrong decoding
    [0,0,1,1,0,1,0,1,1,0,1,0,1,1,0,0,0,1,0,1,0,0,1,1],#Bpt_even ,wrong position wrong decoding
    [0,1,0,1,0,0,1,1,0,0,0,0,0,1,1,0,1,0,0,1,1,1,1,1]]#Bpt_odd  ,correct position wrong distance,correct decoding
    for i in range(4):
        print v[i][0:4],v[i][4:8],v[i][8:12],v[i][12:16],v[i][16:20],v[i][20:24]
        nearestLatticePoint,weight,codeword,winner = decode(leechEncoder(v[i] ))
        print codeword[0:4],codeword[4:8],codeword[8:12],codeword[12:16],codeword[16:20],codeword[20:24]
        print ""
             
        #print '{' ,     
        #for j in range(12):print '{',nearestLatticePoint[j*2],',',nearestLatticePoint[j*2+1],'},',
        #print '}'
    from random import randint
    for i in range(10):
        s = [randint(0,1) for i in range(24)]
        leechEncoder(s)
        print s[0:4],s[4:8],s[8:12],s[12:16],s[16:20],s[20:24]
        nearestLatticePoint,weight,codeword,winner = decode(leechEncoder(s ))
        print codeword[0:4],codeword[4:8],codeword[8:12],codeword[12:16],codeword[16:20],codeword[20:24]
        gcode = multiLevelGolayConstruct(s[0:12])
        print gcode[0:4],gcode[4:8],gcode[8:12],gcode[12:16],gcode[16:20],gcode[20:24]
        print ""
        #print '{' ,     
        #for j in range(12):print '{',nearestLatticePoint[j*2],',',nearestLatticePoint[j*2+1],'},',
        #print '}'
    

def writesomeRandomVectors():
    import random
    l = 2
    errors = []
    
    c = 0
    for t in range(1,1000):
        data = [random.randint(0,1) for i in range(24)] 
        sent= leechEncoder(data)
        newData= multiLevelGolayConstruct(data[0:12])
        k=0
        for d in range(len(newData)):
            k = k+newData[d]*2**d
        
        print data,sent,k

if __name__ == '__main__':
    pass
    #test()
    #test2()
    test3()
    #writesomeRandomVectors()

