/*Author: Lee Carraher
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
#                               1 \/ 3
#                         /\
#                           / 2 \


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
*/

unsigned char xCoords[] = {1,3,5,7,3,5,7,1,3,5,7,1,5,7,1,3};
unsigned char yCoords[] = {7,7,3,3,5,5,1,1,1,1,5,5,7,7,3,3};

inline void printVecI(unsigned char* v,int size){
  unsigned char i = 0;
for(;i<size;i++)printf("%d, ",v[i]);
printf("\n");
}

inline void printVecF(float* v,int size){
unsigned char i = 0;
for(;i<size;i++)printf("%f, ",v[i]);
printf("\n");
}

inline float quicksqrt(float b)
{
    /*
        this thing converges really quickly this is more than enough for fp
    */

    float x = 1.1;
    unsigned char i =0;

    for(;i<25;i++){
        x = (x+(b/x))/2.0;
    }

    return x;
}

inline float distance(float cp[2],float pt[2])
{
    float s = quicksqrt((cp[0]-pt[0])*(cp[0]-pt[0]) + (cp[1]-pt[1])*(cp[1]-pt[1]));
    return s;

}


static void print(unsigned long ret,int ct){
    int i;
    for(i=0;i<ct;i++)
    {
        printf("%d",ret&1);
        ret=ret>>1;
        printf("%d",ret&1);
        ret=ret>>1;
        printf("%d",ret&1);
        ret=ret>>1;
        printf("%d ",ret&1);
        ret=ret>>1;
    } printf("\n");
}


unsigned char H6CodeWords[4][4][4][3]  = {
    {{{0,0,0},{1,1,1},{2,2,2},{3,3,3}},
        {{1,2,3},{0,3,2},{3,0,1},{2,1,0}},
        {{2,3,1},{3,2,0},{0,1,3},{1,0,2}},
        {{3,1,2},{2,0,3},{1,3,0},{0,2,1}}},
    {
        {{1,3,2},{0,2,3},{3,1,0},{2,0,1}},
        {{0,1,1},{1,0,0},{2,3,3},{3,2,2}},
        {{3,0,3},{2,1,2},{1,2,1},{0,3,0}},
        {{2,2,0},{3,3,1},{0,0,2},{1,1,3}}
    },
    {
        {{2,1,3},{3,0,2},{0,3,1},{1,2,0}},
        {{3,3,0},{2,2,1},{1,1,2},{0,0,3}},
        {{0,2,2},{1,3,3},{2,0,0},{3,1,1}},
        {{1,0,1},{0,1,0},{3,2,3},{2,3,2}}
    },
    {
        {{3,2,1},{2,3,0},{1,0,3},{0,1,2}},
        {{2,0,2},{3,1,3},{0,2,0},{1,3,1}},
        {{1,1,0},{0,0,1},{3,3,2},{2,2,3}},
        {{0,3,3},{1,2,2},{2,1,1},{3,0,0}}
    }
};

//000, 110 , 001, 111
float evenAPts[4][2] = {{1.0, 7.0},{5.0, 7.0},{5.0, 3.0},{1.0, 3.0}};
//010 100 011 101
float oddAPts[4][2]  ={{3.0, 5.0},{3.0, 1.0},{7.0, 1.0},{7.0, 5.0}};
//000, 110 , 001, 111
float evenBPts[4][2] = {{3.0, 7.0},{7.0, 7.0},{7.0, 3.0},{3.0, 3.0}};
//010 100 011 101
float oddBPts[4][2]  = {{5.0, 5.0},{5.0, 1.0},{1.0, 1.0},{1.0, 5.0}};

void QAM(float r[12][2], float evenPts[4][2],float oddPts[4][2],float dijs[12][4],float dijks[12][4],unsigned char kparities[12][4]){

//void QAM(float *r, float *evenPts,float *oddPts,float *dijs,float *dijks,int *kparities){
    /*
        this function returns all of the pertinant information from the decoder such as minimum distances, nearest coset leader quadrant, and alternative k-parity distances


    #these maps are seperated into the quadrants of a cartesian plane
    #now we gotta order these properly


    #another simple fix is that the quadrants of QAM be abstractly defined, and the -,+ of order
    #pairs be used to tile the generalized 16bit qam, besides this has to be done anyway so we
    #can get out the real number coordinates in the end
     */

    //the closest even-type Z2 lattice point is used as the
    //coset representatives for all points, not currently used
    //quadrant = [0 for k in range(12)]

    unsigned char i = 0;
    for(;i<12;i++){

        float dist000 = distance(r[i],evenPts[0]);
        float dist110 = distance(r[i],evenPts[1]);
        float dist001 = distance(r[i],evenPts[2]);
        float dist111 = distance(r[i],evenPts[3]);

        if(dist000<dist001)
        {
             dijs[i][0]=dist000;
             dijks[i][0]=dist001;
             kparities[i][0] = 0;
        }
        else{
             dijs[i][0]=dist001;
             dijks[i][0]=dist000;
             kparities[i][0] = 1;
        }
        if(dist110<dist111){
             dijs[i][3]=dist110;
             dijks[i][3]=dist111;
             kparities[i][3] = 0;
        }
        else{
             dijs[i][3]=dist111;
             dijks[i][3]=dist110;
             kparities[i][3] = 1;
        }
        //quadrant[i] = 0


        //min over odds
        float dist010 = distance(r[i],oddPts[0]);
        float dist100 = distance(r[i],oddPts[1]);
        float dist011 = distance(r[i],oddPts[2]);
        float dist101 = distance(r[i],oddPts[3]);
        if (dist010<dist011){
             dijs[i][1]=dist010;
             dijks[i][1]=dist011;
             kparities[i][1] = 0;
        }
        else{
             dijs[i][1]=dist011;
             dijks[i][1]=dist010;
             kparities[i][1] = 1;
        }
        if (dist100<dist101){
             dijs[i][2]=dist100;
             dijks[i][2]=dist101;
             kparities[i][2] = 0;
        }
        else{
             dijs[i][2]=dist101;
             dijks[i][2]=dist100;
             kparities[i][2] = 1;
        }

    }






}



void blockConf(float dijs[12][4],float muEs[6][4],float muOs[6][4],unsigned char prefRepE[6][4][4],unsigned char prefRepO[6][4][4]){
    /*
        computes the Z2 block confidence of the concatonated points projections onto GF4 characters
    */

    //each two symbols is taken as a single character in GF4
    unsigned char i=0;
    for(; i<6;i++){

        //0000 1111
        float s = dijs[2*i][0]+dijs[2*i+1][0];
        float t = dijs[2*i][3]+dijs[2*i+1][3];
        if(s<t){
            muEs[i][0] = s;
            prefRepE[i][0][0] = 0;
            prefRepE[i][0][1] = 0;
            prefRepE[i][0][2] = 0;
            prefRepE[i][0][3] = 0;

        }
        else{
            muEs[i][0] = t;
            //prefRepE[i][0] = 15;//[1,1,1,1]
            prefRepE[i][0][0] = 1;
            prefRepE[i][0][1] = 1;
            prefRepE[i][0][2] = 1;
            prefRepE[i][0][3] = 1;
        }

        //0011 1100 0 3 3 0
        s = dijs[2*i][0]+dijs[2*i+1][3];
        t = dijs[2*i][3]+dijs[2*i+1][0];
        if(s<t){
            muEs[i][1] = s;
            //prefRepE[i][1] = 3;//[0,0,1,1]
            prefRepE[i][1][0] = 0;
            prefRepE[i][1][1] = 0;
            prefRepE[i][1][2] = 1;
            prefRepE[i][1][3] = 1;
        }
        else{
            muEs[i][1] = t;
            //prefRepE[i][1] = 12;//[1,1,0,0]
            prefRepE[i][1][0] = 1;
            prefRepE[i][1][1] = 1;
            prefRepE[i][1][2] = 0;
            prefRepE[i][1][3] = 0;
        }


        //1010 0101
        s = dijs[2*i][2]+dijs[2*i+1][2];
        t = dijs[2*i][1]+dijs[2*i+1][1];
        if (s<t){
            muEs[i][2] = s;
            //prefRepE[i][2] = 10;//[1,0,1,0]
            prefRepE[i][2][0] = 1;
            prefRepE[i][2][1] = 0;
            prefRepE[i][2][2] = 1;
            prefRepE[i][2][3] = 0;
            }
        else{
            muEs[i][2] = t;
            //prefRepE[i][2] = 5;//[0,1,0,1]
            prefRepE[i][2][0] = 0;
            prefRepE[i][2][1] = 1;
            prefRepE[i][2][2] = 0;
            prefRepE[i][2][3] = 1;
        }
        //0110 1001
        s = dijs[2*i][1]+dijs[2*i+1][2];
        t = dijs[2*i][2]+dijs[2*i+1][1];
        if(s<t){
            muEs[i][3] = s;
            //prefRepE[i][3] =6;// [0,1,1,0]
            prefRepE[i][3][0] = 0;
            prefRepE[i][3][1] = 1;
            prefRepE[i][3][2] = 1;
            prefRepE[i][3][3] = 0;
        }
        else{
            muEs[i][3] = t;
            //prefRepE[i][3] = 9;//[1,0,0,1]
            prefRepE[i][3][0] = 1;
            prefRepE[i][3][1] = 0;
            prefRepE[i][3][2] = 0;
            prefRepE[i][3][3] = 1;
        }



    //this operation could be parallel, but probably doesnt need to be

        //1000 0111
        s = dijs[2*i][2]+dijs[2*i+1][0];
        t = dijs[2*i][1]+dijs[2*i+1][3];
        if(s<t){
            muOs[i][0] = s;
            //prefRepO[i][0] = 8;//[1,0,0,0]
            prefRepO[i][0][0] = 1;
            prefRepO[i][0][1] = 0;
            prefRepO[i][0][2] = 0;
            prefRepO[i][0][3] = 0;
        }
        else{
            muOs[i][0] = t;
            //prefRepO[i][0] = 7;//[0,1,1,1]
            prefRepO[i][0][0] = 0;
            prefRepO[i][0][1] = 1;
            prefRepO[i][0][2] = 1;
            prefRepO[i][0][3] = 1;
        }

        //0100 1011
        s = dijs[2*i][1]+dijs[2*i+1][0];
        t = dijs[2*i][2]+dijs[2*i+1][3];
        if (s<t){
            muOs[i][1] = s;
            //prefRepO[i][1] = 4;//[0,1,0,0]
            prefRepO[i][1][0] = 0;
            prefRepO[i][1][1] = 1;
            prefRepO[i][1][2] = 0;
            prefRepO[i][1][3] = 0;
        }
        else{
            muOs[i][1] = t;
            //prefRepO[i][1] = 11;//[1,0,1,1]
            prefRepO[i][1][0] = 1;
            prefRepO[i][1][1] = 0;
            prefRepO[i][1][2] = 1;
            prefRepO[i][1][3] = 1;
        }

        //0010 1101
        s = dijs[2*i][0]+dijs[2*i+1][2];
        t = dijs[2*i][3]+dijs[2*i+1][1];
        if(s<t){
            muOs[i][2] = s;
            //prefRepO[i][2] =2;// [0,0,1,0]
            prefRepO[i][2][0] = 0;
            prefRepO[i][2][1] = 0;
            prefRepO[i][2][2] = 1;
            prefRepO[i][2][3] = 0;
        }
        else{
            muOs[i][2] = t;
            //prefRepO[i][2] = 13;//[1,1,0,1]
            prefRepO[i][2][0] = 1;
            prefRepO[i][2][1] = 1;
            prefRepO[i][2][2] = 0;
            prefRepO[i][2][3] = 1;
        }

        //0001 1110
        s = dijs[2*i][0]+dijs[2*i+1][1];
        t = dijs[2*i][3]+dijs[2*i+1][2];
        if(s<t){
            muOs[i][3] = s;
           // prefRepO[i][3] = 1;//[0,0,0,1]
            prefRepO[i][3][0] = 0;
            prefRepO[i][3][1] = 0;
            prefRepO[i][3][2] = 0;
            prefRepO[i][3][3] = 1;
        }
        else{
            muOs[i][3] = t;
            //prefRepO[i][3] = 14;//[1,1,1,0]
            prefRepO[i][3][0] = 1;
            prefRepO[i][3][1] = 1;
            prefRepO[i][3][2] = 1;
            prefRepO[i][3][3] = 0;
        }
    }

}

void constructHexWord(float mus[6][4],unsigned char chars[6],float charwts[6]){
    /*here we are looking for the least character in the H6 hexacdoe word
       returns the hexacode word and the wt, for using in locating the least reliable symbol
    */
    unsigned char i = 0;
    for(;i<6;i++)
    {
        unsigned char leastChar = 0;
        float leastwt = mus[i][0];

        if(mus[i][1]<leastwt){
            leastwt = mus[i][1];
            leastChar = 1;
        }

        if(mus[i][2]<leastwt){
            leastwt = mus[i][2];
            leastChar = 2;
        }

        if(mus[i][3]<leastwt){
            leastwt = mus[i][3];
            leastChar = 3;
        }

        chars[i] = leastChar;
        charwts[i]=leastwt;

        //printVecI(chars);
        //printVecF(charwts);
    }
}




float minH6(unsigned char  y[6],float charwts[6],float mus[6][4]){
    /*
        this is the minimization over the hexacode function using the 2nd algorithm of  amrani and be'ery ieee may '96
    */


    //locate least reliable
    float leastreliablewt = charwts[0];
    unsigned char leastreliablechar = 0;
    if(charwts[1]>leastreliablewt){
        leastreliablewt = charwts[1];
        leastreliablechar = 1;
    }
    if(charwts[2]>leastreliablewt){
        leastreliablewt = charwts[2];
        leastreliablechar = 2;
    }
    //build candidate list
    unsigned char candslst[8][6]={{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},
                                                        {0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};




    unsigned char i = 0;
    for(;i<4;i++){
        //leastreliable = i

        y[leastreliablechar] = i;

        candslst[i][0] = y[0];
        candslst[i][1] = y[1];
        candslst[i][2] = y[2];
        candslst[i][3] = H6CodeWords[y[0]][y[1]][y[2]][0];

        candslst[i][4] = H6CodeWords[y[0]][y[1]][y[2]][1];
        candslst[i][5] = H6CodeWords[y[0]][y[1]][y[2]][2];
    }

    //y2
    //locate the least reliable symbol in each
    leastreliablewt = charwts[3];
    leastreliablechar = 3;
    if(charwts[4]>leastreliablewt){
        leastreliablewt = charwts[4];
        leastreliablechar = 4;
    }
    if(charwts[5]>leastreliablewt){
        leastreliablewt = charwts[5];
        leastreliablechar = 5;
    }

    i = 0;
    for(;i<4;i++){
        //leastreliable = i
        y[leastreliablechar] = i;
        candslst[i+4][0] = y[3];
        candslst[i+4][1] = y[4];
        candslst[i+4][2] = y[5];
        candslst[i+4][3] = H6CodeWords[y[3]][y[4]][y[5]][0];
        candslst[i+4][4] = H6CodeWords[y[3]][y[4]][y[5]][1];
        candslst[i+4][5] = H6CodeWords[y[3]][y[4]][y[5]][2];
    }

    //minimize over the 8 candidate Hexacode words
    float minCodeWt = 1000.0;
    //minCodeWord = [], this is chars
    i = 0;
    unsigned char j = 0;
    unsigned char  min = 0;
    for(;i<8;i++){
        float m_dist = 0.0;
        j=0;
        for(;j<6;j++)m_dist = m_dist+ mus[j][candslst[i][j]];
        if(m_dist < minCodeWt){
            minCodeWt = m_dist;
            min = i;
        }
    }
    //requires a deep copy here
    for(i=0;i<6;i++)y[i] = candslst[min][i];



    return minCodeWt;
}

float hparity(float weight,unsigned char hexword[6],unsigned char prefReps[6][4][4],float dijs[12][4],unsigned char oddFlag,unsigned char * codeword){
    /*
        here we are resolving the h-parity. which requires that the overall least significant bit parities equal the
        bit parities of each projected GF4 block. aka column parity must equal 1st row parity
    */
   unsigned char parity= 0;
   unsigned char i =0;



    for(;i<6;i++){

        codeword[i*4]=prefReps[i][hexword[i]][0];
        codeword[i*4+1]=prefReps[i][hexword[i]][1];
        codeword[i*4+2]=prefReps[i][hexword[i]][2];
        codeword[i*4+3]=prefReps[i][hexword[i]][3];

        parity = parity + prefReps[i][hexword[i]][0];//this should be the highest order bit
    }

    if((parity&1) == oddFlag)
        return weight;





    float leastwt = 1000.0;
    unsigned char least = 0;
    float deltaX;
    unsigned char idx1,idx2,idxComp1,idxComp2,proj;
    i = 0;
    //walk along the codeword again
    for(;i<6;i++){



        idxComp1 =((codeword[4*i]<<1) +codeword[4*i+1])^3;
        idxComp2 =((codeword[4*i+2]<<1) +codeword[4*i+3])^3;

        deltaX = (dijs[2*i][idxComp1] + dijs[2*i+1][idxComp2]) - (dijs[2*i][idx1] + dijs[2*i+1][idx2]);
        if (deltaX < leastwt){
            leastwt = deltaX;
            least = i*4;
        }
    }


    weight = weight + leastwt;

    codeword[least]= codeword[least]^1;
    codeword[least+1]= codeword[least+1]^1;
    codeword[least+2]= codeword[least+2]^1;
    codeword[least+3]= codeword[least+3]^1;

    return weight;
}

float kparity(float weight,unsigned char * codeword,unsigned char Btype, unsigned char * codeParity, float dijks[12][4], float dijs[12][4],unsigned char kparities[12][4]){
    /*
        this last parity check assures that all A or B points have even/odd parity
    */

    unsigned char parity = 0;
    unsigned char i =0;
    unsigned char idx;
    //int temp = codeword;

    float least =1000;
    float dif;
    unsigned char argLeast = 0;


    for( ;i <12;i++)
     {
        unsigned char n =(codeword[2*i]<<1)+codeword[2*i+1];
         parity= parity+ kparities[i][n];
         codeParity[i] = kparities[i][n];

          dif = dijks[i][n]-dijs[i][n];
          if(dif <= least)
            {
              least = dif;
              argLeast = i;
          }
     }
    if(parity&1 == Btype )
         return weight;

    codeParity[argLeast ]=  codeParity[argLeast ] ^1;
    return weight+least;
}





/*

 //test unencoded
unsigned long ret = 0;

int m= 0 ;
for(;m<12;m++)
  {
    if(r[m][0]>=6)ret = (ret<<2)+3;
    else if(r[m][0]>=4)ret = (ret<<2)+2;
    else if(r[m][0]>=2)ret = (ret<<2)+1;
    else ret = ret<<2;

    if(r[m][1]>=6)ret = (ret<<2)+3;
    else if(r[m][1]>=4)ret = (ret<<2)+2;
    else if(r[m][1]>=2)ret = (ret<<2)+1;
    else ret = ret<<2;

  }
return ret;
 */



//unsigned char* decode(float r[12][2], float *distance){
unsigned long decode(float r[12][2], float *distance){

//

// #####################QAM Dijks ###################
    float dijs[12][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    float dijks[12][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    //there is a set for each quarter decoder, and the A/B_ij odd/even
    unsigned char kparities[12][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    QAM(r,evenAPts,oddAPts,dijs,dijks,kparities);





    // #####################Block Confidences ###################
    //         0  1    w   W
    float muEs[6][4] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}};
    float muOs[6][4] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}};
    unsigned char prefRepE[6][4][4];// = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    unsigned char prefRepO[6][4][4];// = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};

    blockConf(dijs,muEs,muOs,prefRepE,prefRepO);

    unsigned char i;
    /*
    for(i=0;i<12;i++){
      printVecF(r,24);
      printVecF(dijs[i],4);
      printVecF(dijks[i],4);
      printVecI(kparities[i],4);


    }*/



    // #####################Construct Hexacode Word ###################
    unsigned char y[6] = {0,0,0,0,0,0};




    float charwts[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    constructHexWord(muEs,y,charwts);


    // #####################Minimize over the Hexacode ###################
    unsigned char hexword[6] = {0,0,0,0,0,0};
    float weight = minH6(y,charwts,muEs);





    //printf("%d,%d,%d,%d,%d,%d\n",y[0],y[1],y[2],y[3],y[4],y[5]);
    //****chars = y = hexword *****
    unsigned char codeword[24] =  {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
    unsigned char codeParity[12] =  {0,0,0,0, 0,0,0,0, 0,0,0,0} ;





    weight = hparity(weight,y,prefRepE,dijs,0,codeword);//byref
    weight =kparity(weight,codeword,0, codeParity,dijks,dijs,kparities);





    float leastweight = weight;



    unsigned char* leastCodeword = malloc(24*sizeof(unsigned char));
    unsigned long retOpt = 0UL;
    int b;

    for(i=0;i<12;i++){

      //b = (codeword[i*2]<<3) + (codeword[i*2+1]<<2) + (codeParity[i]<<1) + 0;
      //leastCodeword[i*2] = xCoords[b];
      //leastCodeword[i*2+1] = yCoords[b];
      //retOpt = b +(retOpt<<4);

        b = (codeword[i*2]<<2) + (codeword[i*2+1]<<1) + (codeParity[i]);
        retOpt = b +(retOpt<<3);
    }



    //----------------A Odd Quarter Lattice Decoder----------------

    constructHexWord(muOs,y,charwts);
    weight = minH6(y,charwts,muOs);
    weight = hparity(weight,y,prefRepO,dijs,1,codeword);//byref
    weight =kparity(weight,codeword,0,codeParity,dijks,dijs,kparities);


    if(weight<leastweight)
    {
        leastweight = weight;
        retOpt = 0UL;

        for(i=0;i<12;i++){

          //b = (codeword[i*2]<<3) + (codeword[i*2+1]<<2) + (codeParity[i]<<1) + 0;
          //leastCodeword[i*2] = xCoords[b];
          //leastCodeword[i*2+1] = yCoords[b];
          //retOpt = b +(retOpt<<4);

          b = (codeword[i*2]<<2) + (codeword[i*2+1]<<1) + (codeParity[i]);
          retOpt = b +(retOpt<<3);


          /*
          leastCodeword[i*2] = codeword[i*2];
          leastCodeword[i*2+1] = codeword[i*2+1];
          leastCodeword[24+i] = codeParity[i];
          */
        }




    }

    //----------------H_24 Half Lattice Decoder for B points----------------
    QAM(r,evenBPts,oddBPts,dijs,dijks,kparities);
    blockConf(dijs,muEs,muOs,prefRepE,prefRepO);


    //----------------B Even Quarter Lattice Decoder----------------
    constructHexWord(muEs,y,charwts);
    weight = minH6(y,charwts,muEs);

    weight = hparity(weight,y,prefRepE,dijs,0,codeword);//byref
    weight =kparity(weight,codeword,1,codeParity,dijks,dijs,kparities);

    if(weight<leastweight){
        retOpt = 0UL;

        leastweight = weight;




        for(i=0;i<12;i++){

            b = (codeword[i*2]<<2) + (codeword[i*2+1]<<1) + (codeParity[i]);
            retOpt = b +(retOpt<<3);

            //b = (codeword[i*2]<<3) + (codeword[i*2+1]<<2) + (codeParity[i]<<1) + 1;
            //leastCodeword[i*2] = xCoords[b];
            //leastCodeword[i*2+1] = yCoords[b];
            //retOpt = b +(retOpt<<4);
        }

        //for(i=0;i<24;i++)
        //  leastCodeword[i] = codeword[i];
        //for(i=0;i<12;i++)
        //  leastCodeword[i+24] = codeParity[i];

    }

    //----------------B Odd Quarter Lattice Decoder----------------
    constructHexWord(muOs,y,charwts);
    weight = minH6(y,charwts,muOs);
    weight = hparity(weight,y,prefRepO,dijs,1,codeword);//byref
    weight =kparity(weight,codeword,1,codeParity,dijks,dijs,kparities);


    if(weight<leastweight){
        retOpt = 0UL;
        leastweight = weight;

        for(i=0;i<12;i++)
        {
          //b = (codeword[i*2]<<3) + (codeword[i*2+1]<<2) + (codeParity[i]<<1) + 1;
            b = (codeword[i*2]<<2) + (codeword[i*2+1]<<1) + (codeParity[i]);
            retOpt = b +(retOpt<<3);
          //retOpt = b +(retOpt<<4);
          //leastCodeword[i*2] = xCoords[b];
          //leastCodeword[i*2+1] = yCoords[b];
        }

    }


    *distance = leastweight;
    return retOpt;
    //return leastCodeword;


}






float testDist(float r[12][2] , float q[12][2] ){
  float dist = 0;
  int i;
  for (i=0;i<12;i++){
      dist += (r[i][0]- q[i][0]) *(r[i][0]- q[i][0])   +( r[i][1]-  q[i][1])*(r[i][1]-  q[i][1]);
  }

return quicksqrt(dist);
}

void testRatios(int tests , float bias, int * counts, int * totals, float * buckets){
  unsigned int RAND_MAX = 2147483647;

  int i,j;
  float r[12][2];
  float q[12][2];
      for (i=0;i<tests;)
        {

          float dist;
          for(j=0;j<12;j++){
              r[j][0]=6*((float)rand())/RAND_MAX +1;
              r[j][1]=6*((float)rand())/RAND_MAX+1;
          }

          for(j=0;j<12;j++)
            {
              q[j][0]=r[j][0]   +   bias*(-1+2*((float)rand())/RAND_MAX);
              q[j][1]=r[j][1]   +  bias*(-1+2*((float)rand())/RAND_MAX);
          }
          dist = testDist(r,q);

          int b = 0;
          while(dist>buckets[b++]);


          totals[b]++;

          unsigned char * ret1 = decode(r,&dist);

          unsigned char * ret2 = decode(q,&dist);

          if( testDist(r,(float*)ret1)<16.0  &&  testDist(q,(float*)ret2)<16.0)
          {
            unsigned char tr = 1;
            for(j=0;j<24 && tr;j++)
                if (ret1[j]!=ret2[j])tr = 0;
            counts[b]+=tr;

            i++;
          }


      }



}



// The ELFhash function
unsigned long ELFhash(char* key, long M) {

  unsigned long g,h = 0;
  while(*key) {
    h = (h << 4) + *key++;
    g = h & 0xF0000000L;
    if (g) h ^= g >> 24;
    h &= ~g;
  }
  return h % M;
}






void countNZ()
{

  unsigned int RAND_MAX = 2147483647;

   int j,k;
   long i,total=0;

   unsigned char* allhashes = malloc(16777216*32);
   printf("OK\n");

   //zero it out
   for(i=0;i <16777216*32;i++ )allhashes[i]=0;
   printf("OK\n");

   float r[12][2];
   float dist;
   unsigned long ret1;

    for (i=0;i<16777216*16;i++)
    {


           for(j=0;j<12;j++){
               r[j][0]=8*((float)rand())/RAND_MAX;
               r[j][1]=8*((float)rand())/RAND_MAX;
           }

            ret1 = decode(r,&dist);

            ret1 = ELFhash(&ret1, 16777216*32);

           if(allhashes[ret1]  ==0)total++;
           allhashes[ret1]++;
           if(i %100000==0)
           {

               printf("%d:%d\n", i,total);

           }


       }



}





int main(int argc, char* argv[])
{



  srand((unsigned int)12412471);

  countNZ();
  /*
  int i,j,p;
    //some more samples
    //
    //{{7.0,3.0},{3.0,3.0},{7.0,7.0},{3.0,3.0},
    //    {7.0,7.0},{7.0,7.0},{3.0,7.0},{7.0,7.0},
    //    {5.0,5.0},{5.0,1.0},{5.0,5.0},{5.0,5.0}}

    //{{1.0,1.0},{3.0,7.0},{7.0,3.0},{5.0,5.0},
        //{7.0,7.0},{5.0,5.0},{1.0,5.0},{3.0,7.0},
        //{1.0,1.0},{7.0,7.0},{7.0,7.0},{1.0,1.0}}

    //{{5.0,7.0},{3.0,5.0},{1.0,3.0},{7.0,1.0},
        //{3.0,5.0},{5.0,7.0},{3.0,1.0},{5.0,3.0},
        //{1.0,3.0},{3.0,5.0},{1.0,3.0},{7.0,1.0}}

    float dist;
    float r1[12][2] ={{7.0,7.0},{1.0,5.0},{1.0,1.0},{3.0,3.0},
                        {5.0,1.0},{3.0,7.0},{3.0,7.0},{5.0,5.0},
                        {3.0,3.0},{1.0,1.0},{5.0,5.0},{3.0,7.0}};
                    //{{7.0, 5.0}, {3.0, 1.0}, {3.0, 1.0}, {3.0, 1.0}
                    //, {1.0, 3.0}, {1.0, 3.0}, {1.0, 7.0}, {5.0, 3.0},
                    //  {7.0, 5.0}, {7.0, 5.0}, {7.0, 1.0}, {7.0, 1.0}};

    //print(decode(r1,&dist),9);printf("%f\n",dist);
    printVecI(decode(r1,&dist),24);
    long ret = 0;

    //for(i = 0;i<36;i++)ret = decode(r1,&dist)[i] +(ret <<1) ;
    //print(ret,9);


    //print(decode(r1,&dist),9);
    //decode(r1,&dist);
    //add some noise
    r1[0][1]=r1[0][1]-.51;
    r1[3][0]=r1[3][0]+1.01;
    r1[5][1]=r1[5][1]+2.01;
    r1[7][1]=r1[7][1]-.7;
    r1[10][1]=r1[10][1]+1.9;
    r1[11][1]=r1[11][1]+1.9;


    //print(decode(r1,&dist),9);printf("%f\n",dist);
    printVecI(decode(r1,&dist),24);

    //ret = 0;

    //for(i = 0;i<36;i++)ret = decode(r1,&dist)[i] +(ret <<1) ;
    //print(ret,9);

    printf("\n");
    float r2[12][2] =    {{7.0,3.0},{3.0,3.0},{7.0,7.0},{3.0,3.0},
        {7.0,7.0},{7.0,7.0},{3.0,7.0},{7.0,7.0},
        {5.0,5.0},{5.0,1.0},{5.0,5.0},{5.0,5.0}};
                    //{{3.0, 1.0}, {5.0, 7.0}, {7.0, 1.0}, {5.0, 7.0},
                    //  {5.0, 3.0}, {7.0, 5.0}, {1.0, 3.0}, {3.0, 1.0},
                    //  {7.0, 1.0}, {5.0, 3.0}, {7.0, 5.0}, {5.0, 3.0}};

    //print(decode(r2,&dist),9);printf("%f\n",dist);
    printVecI(decode(r2,&dist),24);
    ret = 0;
    //for(i = 0;i<36;i++)ret = (ret <<1) +decode(r2,&dist)[i];
    //print(ret,9);
        //add some noise

        r2[0][1]=r2[0][1]-1.5;
        r2[5][1]=r2[5][1]+2.1;
        r2[6][1]=r2[6][1]-2.1;
        r2[8][1]=r2[8][1]-1.1;
        r2[6][0]=r2[6][0]-2.1;
        r2[8][0]=r2[8][0]-1.1;

      //  print(decode(r2,&dist),9);printf("%f\n",dist);
        printVecI(decode(r2,&dist),24);
        ret = 0;
        //for(i = 0;i<36;i++)ret = (ret <<1) +decode(r2,&dist)[i];
        //print(ret,9);
        printf("\n");


    float r3[12][2] =
                    {{3.0, 3.0}, {3.0, 3.0}, {5.0, 1.0}, {5.0, 5.0},
                    {3.0, 7.0}, {3.0, 3.0}, {5.0, 1.0}, {1.0, 5.0},
                    {3.0, 7.0}, {3.0, 7.0}, {1.0, 5.0}, {5.0, 5.0}};

    //print(decode(r3,&dist),9);printf("%f\n",dist);
    printVecI(decode(r3,&dist),24);
        //add some noise
        r3[0][1]=r3[0][1]-2.5;
        r3[0][0]=0.0;
        r3[3][1]=r3[3][1]+1.51;
        r3[3][0]=r3[3][0]+1.51;
        r3[11][0]=r3[11][0]-1.1;
        r3[11][1]=r3[11][1]-1.1;
        //decode(r3,&dist);
        //print(decode(r3,&dist),9);printf("%f\n",dist);
        printVecI(decode(r3,&dist),24);
        printf("\n");


    float r4[12][2] ={{1.0,1.0},{3.0,7.0},{7.0,3.0},{5.0,5.0},
                      {7.0,7.0},{5.0,5.0},{1.0,5.0},{3.0,7.0},
                      {1.0,1.0},{7.0,7.0},{7.0,7.0},{1.0,1.0}};
                    //{{5.0, 1.0},{3.0, 3.0},{1.0, 1.0},{3.0, 7.0},
                    //{1.0, 1.0},{7.0, 7.0}, {5.0, 5.0}, {3.0, 3.0},
                    //{1.0, 1.0}, {7.0, 3.0}, {1.0, 1.0}, {3.0, 7.0}};

        //print(decode(r4,&dist),9);printf("%f\n",dist);
        printVecI(decode(r4,&dist),24);
        //add some noise
        r4[0][0]=r4[0][0]-1.3;
        r4[1][1]=r4[1][1]+.4;
        r4[2][1]=r4[2][1]-1.3;
        r4[3][0]=r4[3][0]-3.1;
        r4[4][1]=r4[4][1]+2.1;


        r4[11][0]=r4[11][0]-2.1;


        //print(decode(r4,&dist),9);printf("%f\n",dist);
        printVecI(decode(r4,&dist),24);



        //3.9 / 30
        p=30;//this needs a define macro
        float k = 3.9/((float)p);
        float buckets[30] = {0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0 };
        for(i=0;i<p;i++)buckets[i] = ((float)i)*k;


        int counts[30]= {0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0 };
        int totals[30]= {0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0 };





        for(i=1;i<100;i++)
            testRatios(1000 , i*.01 , & counts, & totals,  buckets);


            for (i=0;i<p;i++)printf("%f,",buckets[i]);
            printf("\n");
            for (i=0;i<p;i++)printf("%f,",((float)counts[i])/((float)totals[i]));
            printf("\n");
            for (i=0;i<p;i++)printf("%d,",totals[i]);
            printf("\n");
  */
  return 0;
}


