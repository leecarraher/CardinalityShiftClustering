#Author: Lee Carraher
#Institution: University of Cincinnati, Computer Science Dept.


# this is a nearest lattice point decoder based on the hexacode based decoder of 
#Amrani, Be'ery IEEE Trans. on Comm. '96, with initial construction 
#from  Amrani, Be'ery,Vardy, Sun,Tilborg IEEE Info Thry'94

# the goal is to rewrite this algorithm in efficient C for cuda
# and eventual use as a Hashing Function
# for use in a Cuda Parallel Locality Hash Based Clustering algorthm
# additional implementation may include MPI/Cuda, and 
#anonymous offline data clusering


#-------------QAM Stuff ----------------------
# use a curtailed QAM for all positive signals
#  4 A000 B000 A110 B110
#  3 B101 A010 B010 A101 
#  2 A111 B111 A001 B001 
#  1 B011 A100 B100 A011 
#  0   1    2   3    4
# still gets rotated    \ 4 / 
#		                1 \/ 3
#                         /\ 
#	                    / 2 \	


# leech decoder uses a rotated Z2 lattice, so to find leading cosets
# just find the nearest point in 64QAM, A,B ; odd, even| to the rotated
# input vector
# rotate using the standard 2d rotation transform
#                      [cos x -sin x ]
#                  R = [sin x  cos x ]    cos(pi/4) = sin(pi/4)=1/sqrt(2)
# for faster C implementation use these binary fp constants
# 1/sqrt(2) = cc3b7f669ea0e63f ieee fp little endian
#           = 3fe6a09e667f3bcc ieee fp big endian
#           = 0.7071067811865475244008
#
#v' = v * R
# integer lattice
#
#  4 A000 B000 A110 B110 | A000 B000 A110 B110
#  3 B101 A010 B010 A101 | B101 A010 B010 A101
#  2 A111 B111 A001 B001 | A111 B111 A001 B001
#  1 B011 A100 B100 A011 | B011 A100 B100 A011
#    --------------------|---------------------
# -1 A000 B000 A110 B110 | A000 B000 A110 B110
# -2 B101 A010 B010 A101 | B101 A010 B010 A101
# -3 A111 B111v A001 B001 | A111 B111 A001 B001
# -4 B011 A100 B100 A011 | B011 A100 B100 A011
#even pts {000,110,111,001}
#odd  pts {010,101,100,011}


import math
theta = math.pi/4.0
def mult(t):
    '''
        this rotates the recieved point pi/4 to correspond with the correct tiled qam
    '''
    return  t[0]*math.cos(theta)-t[1]*math.sin(theta),t[0]*math.sin(theta)+t[1]*math.cos(theta)

def distance(cp,pt):
    '''
    simple euclidean distance formula, can be changed if noise is not gaussian
    '''
    return ((cp[0]-pt[0])**2.0 + (cp[1]-pt[1])**2.0)**.5


#access codewords by their information vector indeces H6CodeWords[firstchar][secondchar][thirdchar]
#ex. wW1 - > H6CodeWords[2][3][1] = [0,1,0] => the codeword is wW1010
H6CodeWords = [[[[0,0,0],[1,1,1],[2,2,2],[3,3,3]],[[1,2,3],[0,3,2],[3,0,1],[2,1,0]],[[2,3,1],[3,2,0],[0,1,3],[1,0,2]],[[3,1,2],[2,0,3],[1,3,0],[0,2,1]]],[[[1,3,2],[0,2,3],[3,1,0],[2,0,1]],[[0,1,1],[1,0,0],[2,3,3],[3,2,2]],[[3,0,3],[2,1,2],[1,2,1],[0,3,0]],[[2,2,0],[3,3,1],[0,0,2],[1,1,3]]],[[[2,1,3],[3,0,2],[0,3,1],[1,2,0]],[[3,3,0],[2,2,1],[1,1,2],[0,0,3]],[[0,2,2],[1,3,3],[2,0,0],[3,1,1]],[[1,0,1],[0,1,0],[3,2,3],[2,3,2]]],[[[3,2,1],[2,3,0],[1,0,3],[0,1,2]],[[2,0,2],[3,1,3],[0,2,0],[1,3,1]],[[1,1,0],[0,0,1],[3,3,2],[2,2,3]],[[0,3,3],[1,2,2],[2,1,1],[3,0,0]]]]




#1,2,3,4 quadrants of the Cartesian plane
#probably should be reorganized, see commments at the bottom regarding the construction of the leech lattice vector

#these are rotated and integer centered
#000, 110 , 001, 111
#evenAPts = [[(-1.5,2.0), (-0.5,3.0), (0.5,2.0), (-0.5,1.0)],[ (-3.5,0.0), (-2.5,1.0), (-1.5,0.0), (-2.5,-1.0)],[(-1.5,-2.0), (-0.5,-1.0), (0.5,-2.0), (-0.5,-3.0)],[(0.5,0.0), (1.5,1.0), (2.5,0.0), (1.5,-1.0)]]
#010 100 011 101
#oddAPts  =[[(-0.5,2.0), (0.5,1.0), (1.5,2.0), (0.5,3.0)],[(-2.5,0.0), (-1.5,-1.0), (-0.5,0.0), (-1.5,1.0)],[(-0.5,-2.0), (0.5,-3.0), (1.5,-2.0), (0.5,-1.0)],[(1.5,0.0), (2.5,-1.0), (3.5,0.0), (2.5,1.0)]]
#000, 110 , 001, 111
#evenBPts = [[(-1.0,2.5), (0.0,3.5), (1.0,2.5), (0.0,1.5)],[(-3.0,0.5), (-2.0,1.5), (-1.0,0.5), (-2.0,-0.5)],[(-1.0,-1.5), (0.0,-0.5), (1.0,-1.5), (0.0,-2.5)],[(1.0,0.5), (2.0,1.5), (3.0,0.5), (2.0,-0.5)]]
#010 100 011 101
#oddBPts  = [[(0.0,2.5), (1.0,1.5), (0.0,0.5), (-1.0,1.5)],[(-2.0,0.5), (-1.0,-0.5), (-2.0,-1.5), (-3.0,-0.5)],[(0.0,-1.5), (1.0,-2.5), (0.0,-3.5), (-1.0,-2.5)],[(2.0,0.5), (3.0,-0.5), (2.0,-1.5), (1.0,-0.5)]]


#not rotated
#000, 110 , 001, 111
#evenAPts = [[(.5,3.5),(2.5,3.5),(2.5,1.5),(.5,1.5)],[(-3.5,3.5),(-1.5,3.5),(-1.5,1.5),(-3.5,1.5)],[(-3.5,-.5),(-1.5,-.5),(-1.5,-2.5),(-3.5,-2.5)],[(.5,-.5),(2.5,-.5),(2.5,-2.5),(.5,-2.5)]]
#010 100 011 101
#oddAPts  =[[(1.5,2.5),(1.5,.5),(3.5,.5),(3.5,2.5)],[(-2.5,2.5),(-2.5,.5),(-.5,.5),(-.5,2.5)],[(-2.5,-1.5),(-2.5,-3.5),(-.5,-3.5),(-.5,-1.5)],[(1.5,-1.5),(1.5,-3.5),(3.5,-3.5),(3.5,-1.5)]]
#000, 110 , 001, 111
#evenBPts = [[(1.5,3.5),(3.5,3.5),(3.5,1.5),(1.5,1.5)],[(-2.5,3.5),(-.5,3.5),(-.5,1.5),(-2.5,1.5)],[(-2.5,-.5),(-.5,-.5),(-.5,-2.5),(-2.5,-2.5)],[(1.5,-.5),(3.5,-.5),(3.5,-2.5),(1.5,-2.5)]]
#010 100 011 101

#oddBPts  = [[(2.5,2.5),(2.5,.5),(.5,.5),(.5,2.5)],[(-1.5,2.5),(-1.5,.5),(-3.5,.5),(-3.5,2.5)],[(-1.5,-1.5),(-1.5,-3.5),(-3.5,-3.5),(-3.5,-1.5)],[(2.5,-1.5),(2.5,-3.5),(.5,-3.5),(.5,-1.5)]]

#for i in range(4):
#    evenAPts[i]=[mult(p) for p in evenAPts[i]]
#    oddAPts[i]=[mult(p) for p in oddAPts[i]]
#    evenBPts[i]=[mult(p) for p in evenBPts[i]]
#    oddBPts[i]=[mult(p) for p in oddBPts[i]]


#x2 all pts
#000, 110 , 001, 111
evenAPts = [[(1.0, 7.0),(5.0, 7.0),(5.0, 3.0),(1.0, 3.0)],[(-7.0, 7.0),(-3.0, 7.0),(-3.0, 3.0),(-7.0, 3.0) ],[(-7.0, -1.0),(-3.0, -1.0),(-3.0, -5.0),(-7.0, -5.0)],[(1.0, -1.0),(5.0, -1.0),(5.0, -5.0),(1.0, -5.0) ]]
#010 100 011 101
oddAPts  =[[(3.0, 5.0),(3.0, 1.0),(7.0, 1.0),(7.0, 5.0)],[(-5.0, 5.0),(-5.0, 1.0),(-1.0, 1.0),(-1.0, 5.0)],[(-5.0, -3.0),(-5.0, -7.0),(-1.0, -7.0),(-1.0, -3.0) ],[(3.0, -3.0),(3.0, -7.0),(7.0, -7.0),(7.0, -3.0) ]]
#000, 110 , 001, 111
evenBPts = [[(3.0, 7.0),(7.0, 7.0),(7.0, 3.0),(3.0, 3.0) ],[(-5.0, 7.0),(-1.0, 7.0),(-1.0, 3.0),(-5.0, 3.0) ],[(-5.0, -1.0),(-1.0, -1.0),(-1.0, -5.0),(-5.0, -5.0) ],[(3.0, -1.0),(7.0, -1.0),(7.0, -5.0),(3.0, -5.0)]]
#010 100 011 101
oddBPts  = [[(5.0, 5.0),(5.0, 1.0),(1.0, 1.0),(1.0, 5.0)],[(-3.0, 5.0),(-3.0, 1.0),(-7.0, 1.0),(-7.0, 5.0) ],[(-3.0, -3.0),(-3.0, -7.0),(-7.0, -7.0),(-7.0, -3.0) ],[(5.0, -3.0),(5.0, -7.0),(1.0, -7.0),(1.0, -3.0) ]]

def simpleQAM(r):
    ret = []
    for i in range(12):
        least = 1000.0
        for k in range(4):
            for j in range(4):
                distance(r[i],evenAPts[k][j])
                if d < least:
                    least = d
                    lret = evenAPts[k][j]
            for j in range(4):
                distance(r[i],oddAPts[k][j])
                if d < least:
                    least = d
                    lret = oddAPts[k][j]
            for j in range(4):
                distance(r[i],evenBPts[k][j])
                if d < least:
                    least = d
                    lret = evenBPts[k][j]
            for j in range(4):
                distance(r[i],oddBPts[k][j])
                if d < least:
                    least = d
                    lret = oddBPts[k][j]
        ret.append(lret)
    return ret      


def QAM(r,evenPts,oddPts):
    '''
        this function returns all of the pertinant information from the decoder such as minimum distances, nearest coset leader quadrant, and alternative k-parity distances
    '''

    #these maps are seperated into the quadrants of a cartesian plane
    #now we gotta order these properly
    #this is a better way to organize these !!!! 00 01 10 11 
    #because of easy parity-ing, too late for the python version, will be fixed in C/cuda
    
    #another simple fix is that the quadrants of QAM be abstractly defined, and the -,+ of order
    #pairs be used to tile the generalized 16bit qam, besides this has to be done anyway so we
    #can get out the real number coordinates in the end
     

    #r is length 12 from R^2 points
    #        00 11 01 10
    dijs = [[[0,0],[0,0]] for k in range(12)]
    dijks = [[[0,0],[0,0]] for k in range(12)]
    #there is a set for each half decoder, and the A/B_ij point selected
    kparity = [[[0,0],[0,0]] for k in range(12)]

    #the closest even-type Z2 lattice point is used as the 
    #coset representatives for all points
    quadrant = [0 for k in range(12)]


    for i in range(12):
                
        dist000 = distance(r[i],evenPts[0][0])
        dist110 = distance(r[i],evenPts[0][1])
        dist001 = distance(r[i],evenPts[0][2])
        dist111 = distance(r[i],evenPts[0][3])
             
        
        if dist000<dist001:
             dijs[i][0][0]=dist000
             dijks[i][0][0]=dist001
             kparity[i][0][0] = 0
        else:
             dijs[i][0][0]=dist001
             dijks[i][0][0]=dist000
             kparity[i][0][0] = 1

        if dist110<dist111:
             dijs[i][0][1]=dist110
             dijks[i][0][1]=dist111
             kparity[i][0][1] = 0
        else:
             dijs[i][0][1]=dist111
             dijks[i][0][1]=dist110
             kparity[i][0][1] = 1
        quadrant[i] = 0
        '''
        for j in xrange(1,4):
            dist000 = distance(r[i],evenPts[j][0])
            dist110 = distance(r[i],evenPts[j][1])
            dist001 = distance(r[i],evenPts[j][2])
            dist111 = distance(r[i],evenPts[j][3])

            if dist000<dijs[i][0][0] or dist001<dijs[i][0][0]:
            
                if dist000<dist001:
                     dijs[i][0][0]=dist000
                     dijks[i][0][0]=dist001
                     kparity[i][0][0] = 0
                else:
                     dijs[i][0][0]=dist001
                     dijks[i][0][0]=dist000
                     kparity[i][0][0] = 1
                quadrant[i] = j

            if dist110<dijs[i][0][1] or dist111<dijs[i][0][1]:
                if dist110<dist111:
                     dijs[i][0][1]=dist110
                     dijks[i][0][1]=dist111
                     kparity[i][0][1] = 0
                else:
                     dijs[i][0][1]=dist111
                     dijks[i][0][1]=dist110
                     kparity[i][0][1] = 1

        '''
        #min over odds
        dist010 = distance(r[i],oddPts[0][0])
        dist100 = distance(r[i],oddPts[0][1])
        dist011 = distance(r[i],oddPts[0][2])
        dist101 = distance(r[i],oddPts[0][3])

        if dist010<dist011:
             dijs[i][1][0]=dist010
             dijks[i][1][0]=dist011
             kparity[i][1][0] = 0
        else:
             dijs[i][1][0]=dist011
             dijks[i][1][0]=dist010
             kparity[i][1][0] = 1      

        if dist100<dist101:
             dijs[i][1][1]=dist100
             dijks[i][1][1]=dist101
             kparity[i][1][1] = 0
        else:
             dijs[i][1][1]=dist101
             dijks[i][1][1]=dist100
             kparity[i][1][1] = 1

    
        '''
        for j in xrange(1,4):
            dist000 = distance(r[i],evenPts[j][0])
            dist110 = distance(r[i],evenPts[j][1])
            dist001 = distance(r[i],evenPts[j][2])
            dist111 = distance(r[i],evenPts[j][3])

            if dist010<dijs[i][1][0] or dist011<dijs[i][1][0]:

                if dist010<dist011:
                     dijs[i][1][0]=dist010
                     dijks[i][1][0]=dist011
                     kparity[i][1][0] = 0
                else:
                     dijs[i][1][0]=dist011
                     dijks[i][1][0]=dist010
                     kparity[i][1][0] = 1

            if dist100<dijs[i][1][1] or dist101<dijs[i][1][1]:
                if dist100<dist101:
                     dijs[i][1][1]=dist100
                     dijks[i][1][1]=dist101
                     kparity[i][1][1] = 0
                else:
                     dijs[i][1][1]=dist101
                     dijks[i][1][1]=dist100
                     kparity[i][1][1] = 1

        '''
        #i think replacement work'd - all codewords are of weight 8,12,16, i guess 0 or 24 too, but ihavent seen any yet
        #Adijs [12 sets]: each contain 4 distances [00,11],[01,10]
        #Adijks [12 sets]: each contain 4 distances [00,11],[01,10] with complement K-parity bits
        #for j in xrange(1,4):
        #    dist010 = distance(r[i],oddPts[j][0])
        #    dist100 = distance(r[i],oddPts[j][1])
        #    dist011 = distance(r[i],oddPts[j][2])
        #    dist101 = distance(r[i],oddPts[j][3])
        #odds
        #    if dijs[i][1][0]>min(dist010,dist011):
        #        dijs[i][1][0]=min(dist010,dist011)
        #        dijks[i][1][0]=max(dist010,dist011)

        #    if dijs[i][1][1]>min(dist100,dist101):
        #        dijs[i][1][1]=min(dist100,dist101) 
        #        dijks[i][1][1]=max(dist010,dist011)
    return dijs,dijks,kparity,quadrant

#dijs is length 12, [[00 11],[ 01 10]]
def blockConf(dijs):
    '''
        computes the Z2 block confidence of the concatonated points projections onto GF4 characters
    '''
    #         0  1    w   W
    muEs = [[0.0,0.0,0.0,0.0] for i in range(6)]
    muOs = [[0.0,0.0,0.0,0.0] for i in range(6)]
    #store the preferable binary rep of the hexacodewords, this should be made into a bitvector
    prefRepE = [[[0,0,0,0]for i in range(4)] for i in range(6)]
    prefRepO = [[[0,0,0,0]for i in range(4)] for i in range(6)]
    #each two symbols is taken as a single character in GF4
    for i in xrange(6):
        
        #0000 1111 00 00 01 01
        s = dijs[2*i][0][0]+dijs[2*i+1][0][0]
        t = dijs[2*i][0][1]+dijs[2*i+1][0][1]
        if s<t:
            muEs[i][0] = s
            prefRepE[i][0] = [0,0,0,0]
        else:
            muEs[i][0] = t
            prefRepE[i][0] = [1,1,1,1]
        
        #0011 1100 00 01 01 00
        s = dijs[2*i][0][0]+dijs[2*i+1][0][1]
        t = dijs[2*i][0][1]+dijs[2*i+1][0][0]
        if s<t:
            muEs[i][1] = s
            prefRepE[i][1] = [0,0,1,1]
        else:
            muEs[i][1] = t
            prefRepE[i][1] = [1,1,0,0]


        #1010 0101 11 11 10 10
        s = dijs[2*i][1][1]+dijs[2*i+1][1][1]
        t = dijs[2*i][1][0]+dijs[2*i+1][1][0]
        if s<t:
            muEs[i][2] = s
            prefRepE[i][2] = [1,0,1,0]
        else:
            muEs[i][2] = t
            prefRepE[i][2] = [0,1,0,1]

        #0110 1001 10 11 11 10
        s = dijs[2*i][1][0]+dijs[2*i+1][1][1]
        t = dijs[2*i][1][1]+dijs[2*i+1][1][0]
        if s<t:
            muEs[i][3] = s
            prefRepE[i][3] = [0,1,1,0]
        else:
            muEs[i][3] = t
            prefRepE[i][3] = [1,0,0,1]



    #this operation could be parallel, but probably doesnt need to be

        #1000 0111 11 00 10 01
        s = dijs[2*i][1][1]+dijs[2*i+1][0][0]
        t = dijs[2*i][1][0]+dijs[2*i+1][0][1]
        if s<t:
            muOs[i][0] = s
            prefRepO[i][0] = [1,0,0,0]
        else:
            muOs[i][0] = t
            prefRepO[i][0] = [0,1,1,1]


        #0100 1011 10 00 11 01
        s = dijs[2*i][1][0]+dijs[2*i+1][0][0]
        t = dijs[2*i][1][1]+dijs[2*i+1][0][1]
        if s<t:
            muOs[i][1] = s
            prefRepO[i][1] = [0,1,0,0]
        else:
            muOs[i][1] = t
            prefRepO[i][1] = [1,0,1,1]

        #0010 1101 00 11 01 10
        s = dijs[2*i][0][0]+dijs[2*i+1][1][1]
        t = dijs[2*i][0][1]+dijs[2*i+1][1][0]
        if s<t:
            muOs[i][2] = s
            prefRepO[i][2] = [0,0,1,0]
        else:
            muOs[i][2] = t
            prefRepO[i][2] = [1,1,0,1]

        #0001 1110 00 10 01 11
        s = dijs[2*i][0][0]+dijs[2*i+1][1][0]
        t = dijs[2*i][0][1]+dijs[2*i+1][1][1]
        if s<t:
            muOs[i][3] = s
            prefRepO[i][3] = [0,0,0,1]
        else:
            muOs[i][3] = t
            prefRepO[i][3] = [1,1,1,0]
        #debug looks like block confidences are ok        
        #print 'evens'        
        #print muEs[i]
        #print 'odds'
        #print muOs[i]
    return muEs,muOs,prefRepE,prefRepO


def constructHexWord(mus):
    '''here we are looking for the least character in the H6 hexacdoe word
       returns the hexacode word and the wt, for using in locating the least reliable symbol
    '''
    chars = [0,0,0,0,0,0]
    charwts = [0.0,0.0,0.0,0.0,0.0,0.0]
    for i in xrange(6):
        leastChar = 0
        leastwt = mus[i][0]

        if mus[i][1]<leastwt:
            leastwt = mus[i][1]
            leastChar = 1

        if mus[i][2]<leastwt:
            leastwt = mus[i][2]
            leastChar = 2

        if mus[i][3]<leastwt:
            leastwt = mus[i][3]
            leastChar = 3

        
        chars[i] = leastChar
        charwts[i]=leastwt
    return chars,charwts


#
#
def minH6(chars,charwts,mus):
    '''
        this is the minimization over the hexacode funtion using the 2nd algorithm of  amrani and be'ery ieee may '96
    '''
    #split the codeword into y1 & y2
    #y1
    #locate least reliable
    leastreliablewt = charwts[0]
    leastreliablechar = 0
    if charwts[1]>leastreliablewt:
        leastreliablewt = charwts[1]
        leastreliablechar = 1
    if charwts[2]>leastreliablewt:
        leastreliablewt = charwts[2]
        leastreliablechar = 2

    #build candidate list
    candslst=[]
    y = chars[0:3]
    y[leastreliablechar] = 0
    candslst.append(y+H6CodeWords[y[0]][y[1]][y[2]])
    y[leastreliablechar] = 1
    candslst.append(y+H6CodeWords[y[0]][y[1]][y[2]])
    y[leastreliablechar] = 2
    candslst.append(y+H6CodeWords[y[0]][y[1]][y[2]])
    y[leastreliablechar] = 3
    candslst.append(y+H6CodeWords[y[0]][y[1]][y[2]])
    

    #y2
    #locate the least reliable symbol in each
    leastreliablewt = charwts[3]
    leastreliablechar = 0
    if charwts[4]>leastreliablewt:
        leastreliablewt = charwts[4]
        leastreliablechar = 1
    if charwts[5]>leastreliablewt:
        leastreliablewt = charwts[5]
        leastreliablechar = 2
    

    #still a little concerned that this should return the same codewords if the 
    #error is contained in the first part one way to do this is to reverse the lookup
    #for the second part of the H6 ie 311 -> 313, then 313 -> 311 for the second half
    #not sure if this is how it is supposed to work, could always check overall error rate at the end
    #with both options subbed in. one way to get this is add y to the beginning and revers idx from 
    #H6[y_0 y_1 y_2] to  rev(H6[y_2 y_1 y_0])
    '''causes 83 errors over 79, is this correct or just random

    y = chars[3:6]
    y[leastreliablechar] = 0
    candslst.append(y+[i for i in reversed(H6CodeWords[y[2]][y[1]][y[0]])   ])
    y[leastreliablechar] = 1
    candslst.append(y+[i for i in reversed(H6CodeWords[y[2]][y[1]][y[0]])  ])
    y[leastreliablechar] = 2
    candslst.append(y+[i for i in reversed(H6CodeWords[y[2]][y[1]][y[0]])  ])
    y[leastreliablechar] = 3
    candslst.append(y+[i for i in reversed(H6CodeWords[y[2]][y[1]][y[0]])  ])
    '''

    
    y = chars[3:6]
    y[leastreliablechar] = 0
    candslst.append(H6CodeWords[y[0]][y[1]][y[2]]+y)
    y[leastreliablechar] = 1
    candslst.append(H6CodeWords[y[0]][y[1]][y[2]]+y)
    y[leastreliablechar] = 2
    candslst.append(H6CodeWords[y[0]][y[1]][y[2]]+y)
    y[leastreliablechar] = 3
    candslst.append(H6CodeWords[y[0]][y[1]][y[2]]+y)
    #'''



    #minimize over the 8 candidate Hexacode words
    minCodeWt = 1000.0
    minCodeWord = []
    #print candslst

    #kmin = 0
    for k in xrange(8):
        m_dist = 0.0
        for j in xrange(6):
            m_dist = m_dist+ mus[j][candslst[k][j]]
        if m_dist < minCodeWt:
            #kmin = k 
            minCodeWt = m_dist
            minCodeWord = candslst[k]

    return minCodeWt, minCodeWord



#00->[0][0]  ~ [0][1]
#01->[1][0]  ~ [1][1]
#10->[1][1]  ~ [1][0]
#11->[0][1]  ~ [0][0]
normal = {(0,0):(0,0),(0,1):(1,0),(1,0):(1,1),(1,1):(0,1)}
complement = {(0,0):(0,1),(0,1):(1,1),(1,0):(1,0),(1,1):(0,0)}

#problem with codeword is not here
def hparity(weight,hexword,prefReps,dijs,oddFlag):
    '''
        here we are resolving the h-parity. which requres that the overall least significant bit parities equal the 
        bit parities of each projected GF4 block. aka column parity must equal 1st row parity
    '''
    
    #all references to codeword are for testing
    codeword = []

    parity = 0
    for i in xrange(6):
        parity = parity + prefReps[i][hexword[i]][0]
        codeword.extend(prefReps[i][hexword[i]])

    #print codeword, parity
    parity = parity%2

    if parity == oddFlag:
        #print "good paritys", weight, codeword
        #print weight
        return weight,codeword

    #print hexword
    #find min parity swap that will fix parity mismatch
    
	#for 6 chars
	#find min{ (dij[0]+dij[1]) - (~dij[0]+~dij[1])}
	#replace with ~dij[0]+~dij[1]
	#weight = weight +[(dij[0]+dij[1]) - (~dij[0]+~dij[1])]
    
    leastwt = 1000.0
    least = 0
    for i in xrange(6):
        proj = prefReps[0][hexword[0]]#this could also be acquired from the above codeword
        #in C version lets make sure to make these correctly correspond so we can just grab them, for now map!
        #ie 00->[0][0],01->[0][1] and ~00->[1][1], ~10->[0][1]
        
        projN = normal[(proj[0],proj[1])]+normal[(proj[2],proj[3])]#this makes a list
        #print proj
        #print projN
        projC = complement[(proj[0],proj[1])]+complement[(proj[2],proj[3])]
        deltaX = (dijs[2*i][projC[0]][projC[1]] + dijs[2*i+1][projC[2]][projC[3]]) - (dijs[2*i][projN[0]][projN[1]] + dijs[2*i+1][projN[2]][projN[3]])
        if deltaX < leastwt:
            leastwt = deltaX
            least = i

    ##############Testing Start################
    #tot = 0
    #for i in xrange(24):tot = tot+codeword[i]
    #if tot!=8 and tot!=16 and tot!=12 and tot!=0 and tot!=24:print codeword, hexword
    #import copy
    #prevCode= copy.deepcopy(codeword)
    #wt =tot
    #print 'after'
    ##############Testing End################


    ### bad ass bitwise way to do this : c ^ (0xF<<((l-1)<<2))
    #update the codeword and complement it
    weight = weight + leastwt
    for i in range(4*least,4*least+4):
        codeword[i] = codeword[i]^1
            

    ##############Testing Start################
    #tot = 0
    #for i in xrange(24):tot = tot+codeword[i]
    #if tot!=8 and tot!=16 and tot!=12 and tot!=0 and tot!=24:
    #    print hexword, range(least-1,least+3)
    #    print '\t'+str(codeword)+' '+str(tot)
    #    print '\t'+str(prevCode)+' '+str(wt)
    ##############Testing End################


    return weight, codeword


def kparity(weight,codeword,dijks,dijs,kparities,Btype):
    '''
        this last parity check assures that all A or B points have even/odd parity
    '''
    #r is length 12 from R^2 points
    #        00 11 01 10
    #dijs = [[[0,0],[0,0]] for k in range(12)]
    #dijks = [[[0,0],[0,0]] for k in range(12)]
    #there is a set for each half decoder, and the A/B_ij point selected
    kparity = [[[0,0],[0,0]] for k in range(12)]
    outparities = [0,0,0,0,0,0,0,0,0,0,0,0]


    parity = 0
    for i in range(12):
        proj = normal[(codeword[i*2],codeword[i*2+1])]
        outparities[i] = kparities[i][proj[0]][proj[1]]
        parity= parity+ outparities[i]
    
    #if btype k-parity is odd
    if parity%2 ==Btype :
        return weight,codeword,outparities


    least = 1000
    leastsymbol = 0
    for i in range(12):
        proj = normal[(codeword[i*2],codeword[i*2+1])]
        dif = dijs[i][proj[0]][proj[1]] - dijks[i][proj[0]][proj[1]]
        if dif < least:
            least = dif
            leastsymbol = i

    outparities[leastsymbol] = outparities[leastsymbol]^1
    weight = weight+least
    return weight,codeword,outparities
        




from ctypes import *
import os


leechdec = cdll.LoadLibrary(os.getcwd()+'/leechDecoder.so')
leechdec.decode.argtypes= [POINTER(c_float),POINTER(c_float)]
leechdec.decode.restype = c_ulong
floatArray =c_float*24

def dec(pts):

    scale(pts,0.,1.,0.0,8.0)
    cx=floatArray(*pts)
    d = c_float()
    result = leechdec.decode(cx,byref(d))
    #r = long(result) & 0xffffff
    #if not (cnt(r)==8 or cnt(r)==12 or cnt(r)==16 or cnt(r)==24) :  print long(result)
    return long(result)




def flatten(l):
    out = []
    for s in l:
        out.extend(s)
    return out

def QAMC(r):

    #stuff for testing C-code
    dijsC = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];
    dijksC = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];
        #there is a set for each quarter decoder, and the A/B_ij odd/even
    kparitiesC = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];
    #x2 all pts
    #000, 110 , 001, 111
    evenAPtsC = [[1.0, 7.0],[5.0, 7.0],[5.0, 3.0],[1.0, 3.0]]
    #010 100 011 101
    oddAPtsC  = [[3.0, 5.0],[3.0, 1.0],[7.0, 1.0],[7.0, 5.0]]
    #000, 110 , 001, 111
    evenBPtsC = [[3.0, 7.0],[7.0, 7.0],[7.0, 3.0],[3.0, 3.0] ]
    #010 100 011 101
    oddBPtsC  = [[5.0, 5.0],[5.0, 1.0],[1.0, 1.0],[1.0, 5.0]]

    r = flatten(r)  
    floatArray = c_float*24
    r=floatArray(*r)
    
    evenAPtsC = flatten(evenAPtsC)
    floatArray =c_float*8
    evenAPtsC = floatArray(*evenAPtsC)
    
    oddAPtsC = flatten(oddAPtsC)
    floatArray =c_float*8
    oddAPtsC = floatArray(*oddAPtsC)

    #dijsC = flatten(dijsC)
    floatArray =(c_float*48)()
    dijsC = cast( floatArray,POINTER(c_float))
    
    #dijksC = flatten(dijksC)
    floatArray =(c_float*48)()
    dijksC = cast(floatArray,POINTER(c_float))

    kparitiesC = flatten(kparitiesC)
    intArray =(c_int*48)()
    kparitiesC = cast(intArray,POINTER(c_int)) #pointer(kparitiesC)
    

    leechdec.QAM(r, evenAPtsC,oddAPtsC,dijsC,dijksC,kparitiesC)
    




def decode(R):
    '''
        this is the decoder for the 4 quarter lattice decoders of even and odd points and A and B points
    '''

    #----------------H_24 Half Lattice Decoder for A points----------------
    dijs,dijks,kparities,quadrant =  QAM(R,evenAPts,oddAPts)
    
     
    #float dijs[12][4],float dijks[12][4],int kparities[12][4]


    muEs,muOs,prefE,prefO = blockConf(dijs)

    #avg( sum([sum(i) for i in muEs])/sum([sum(i) for i in muOs]) )= 1.0
    #this strongly suggest that QAM is returning proper distances

    #----------------A Even Quarter Lattice Decoder----------------
    chars,charwts= constructHexWord(muEs)
    weight ,hexword= minH6(chars,charwts,muEs)
    weight,codeword = hparity(weight,hexword,prefE,dijs,False)
    leastweight,outputCodeword,leastoutparities = kparity(weight,codeword,dijks,dijs,kparities,False)
    Apoint = True
    winner = 0

    #----------------A Odd Quarter Lattice Decoder----------------
    chars,charwts= constructHexWord(muOs)
    weight ,hexword=  minH6(chars,charwts,muOs) 
    weight,codeword = hparity(weight,hexword,prefO,dijs,True)
    weight,codeword,outparities = kparity(weight,codeword,dijks,dijs,kparities,False)
    if weight < leastweight:
        leastweight = weight
        outputCodeword = codeword
        leastoutparities = outparities
        winner = 1
    
    #----------------H_24 Half Lattice Decoder for B points----------------
    dijs,dijks,kparities,quadrant =  QAM(R,evenBPts,oddBPts)
    muEs,muOs,prefE,prefO = blockConf(dijs)


    #----------------B Even Quarter Lattice Decoder----------------
    chars,charwts= constructHexWord(muEs)
    weight ,hexword= minH6(chars,charwts,muEs)
    weight,codeword = hparity(weight,hexword,prefE,dijs,False)
    weight,codeword,outparities = kparity(weight,codeword,dijks,dijs,kparities,True)

    if weight < leastweight:
        outputCodeword = codeword
        leastweight = weight
        leastoutparities = outparities
        Apoint = False
        winner = 2

    #----------------B Odd Quarter Lattice Decoder----------------
    chars,charwts= constructHexWord(muOs)
    weight ,hexword=  minH6(chars,charwts,muOs)
    weight,codeword = hparity(weight,hexword,prefO,dijs,True)
    weight,codeword,outparities = kparity(weight,codeword,dijks,dijs,kparities,True)

    if weight < leastweight:
        outputCodeword = codeword
        leastweight = weight
        leastoutparities = outparities
        Apoint = False
        winner = 3

    ###Something is wrong with the parity selection of the GF4 representations FIXED

    nearestLatticePoint = []
    #lets create a better point set organized by 3 bit binary rep for the C version
    ptSet = [evenAPts,oddAPts]
    if not Apoint: 
        ptSet = [evenBPts,oddBPts]


    #i admit this is totally wonky, the lattice pts vectors need to be re-organized, but for now this works
    for j in range(12):
        # first lets figure out if we are dealing with an even or odd point      
        pt = ptSet[(outputCodeword[j*2]+outputCodeword[j*2+1])%2]
        #now the quadrant
        pt = pt[quadrant[j]]
        #next find the parity part
        index = leastoutparities[j]*2
        #finally which one is it
        index = outputCodeword[j*2]+index

        nearestLatticePoint.extend(pt[index])


    return nearestLatticePoint,leastweight,outputCodeword,winner


#parity correct test! we need more for full proof, but not for use in a LSH algorithm
def simpleQAM(rec):
    ret = []
    for pt in rec:
        leastdist = 1000.0
        leastpt = (0,0)
        for i in range(4):
            for p in evenAPts[0]:
                if  distance(pt,p) < leastdist:
                    leastdist=distance(pt,p)
                    leastpt = p

            for p in oddAPts[0]:
                if  distance(pt,p) < leastdist:
                    leastdist=distance(pt,p)
                    leastpt = p
            for p in evenBPts[0]:
                if  distance(pt,p) < leastdist:
                    leastdist=distance(pt,p)
                    leastpt = p
            for p in oddBPts[0]:
                if  distance(pt,p) < leastdist:
                    leastdist=distance(pt,p)
                    leastpt = p
        ret.append(leastpt[0])
        ret.append(leastpt[1])
    return ret



def simpleTest():
    import random
    import leechencoder
    dataBitVec = random.randint(0,2**24-1)
    data = [int(i) for i in bin(dataBitVec).split('b')[1]]
    while len(data)<24:data = [0]+data

    #this is a lattice point center
    sent= leechencoder.leechEncoder(data)
    
    
    latpts,weight,initcodeword,winningdecoder = decode(sent)
    print latpts
    print initcodeword[0:4],initcodeword[4:8],initcodeword[8:12],initcodeword[12:16],initcodeword[16:20],initcodeword[20:24]
    initcodeword = [int(i) for i in bin(dec(flatten(sent))).split("b")[1]]
    l = 24-len(initcodeword)
    initcodeword = [0]*l + initcodeword
    print initcodeword[20:24],initcodeword[16:20],initcodeword[12:16],initcodeword[8:12],initcodeword[4:8],initcodeword[0:4]
    #altCw = [int(i) for i in bin(dec(flatten(received))).split('b')[1]]



def testLatticeDecoder():

    '''
    import leechencoder
    gwords = []
    for i in range(0,4096):
        strWord = bin(i).split('b')[1]
        l = len(strWord)
        s = [0 for i in range(12-l)]
        for j in range(l):
            s.append(int(strWord[j]))
        gwords.append(leechencoder.multiLevelGolayConstruct(s))


    print "words generated"
    '''
    import random
    import leechencoder
    l = 5000
    errors = []
    errors2 = []
    weights = [0 for i in range(25)]
    decErrors = [0,0,0,0]
    
    for e in range(2,60):
        c = 0
        c2 = 0
        
        for t in range(0,l):
            
            toterr = 0
            data = [random.randint(0,1) for i in range(24)] 
            sent= leechencoder.leechEncoder(data)
            
            
            
            received =[]
            
            pointErrors = random.randint(1,12)
            overallDistance = (float(e)/10.0)/float(pointErrors)
            errlst = []

            for i in range(12):
                received.append((sent[i][0],sent[i][1]))
            for i in range(pointErrors):
                pl = random.randint(0,11)
                while errlst.__contains__(pl):
                    pl = random.randint(0,11)

                errlst.append(pl)
                received[pl] = (received[pl][0]+overallDistance/2.0,received[pl][1]-overallDistance/2.0)


            # some random errors
            
            received[3] = (received[3][0]+overallDistance*24,received[3][1])
            #    received.append((sent[i][0]+overallDistance,sent[i][1]+overallDistance))
        
            latpts,weight,codeword,winningDecoder = decode(received)
            
            weights[sum(codeword)]=weights[sum(codeword)]+1
            newData= leechencoder.multiLevelGolayConstruct(data[0:12])

            if sum([(latpts[2*i]-sent[i][0])+(latpts[2*i+1]-sent[i][1]) for i in range(12)]):
                c = c+1
                decErrors[winningDecoder] = decErrors[winningDecoder]+1
            #if sum([newData[i]-codeword[i] for i in range(24)]):c2 = c2+1

            latpts=simpleQAM(received)
            if sum([(latpts[2*i]-sent[i][0])+(latpts[2*i+1]-sent[i][1]) for i in range(12)]):c2 = c2+1
            


        #if c!=0:print e,overallDistance*pointErrors,c,c2
        #print e
        
        errors.append(c/float(l))
        errors2.append(c2/float(l))
    print decErrors
    import pylab
    import math
    pylab.plot([(float(e)/10.0) for e in range(len(errors))],[math.log(float(errors2[e]/float(errors[e]))) for e in range(len(errors))  ])
    #pylab.plot([(float(e)/10.0) for e in range(len(errors2))],[1.0-e for e in errors2])
    #pylab.legend(['Aeven','Aodd','Beven','Bodd'])
    pylab.show()
    



if __name__ == '__main__':
    pass
    testLatticeDecoder()
    simpleTest()





