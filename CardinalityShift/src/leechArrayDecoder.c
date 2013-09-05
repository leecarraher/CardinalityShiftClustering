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
#  7 A000 B000 A110 B110
#  5 B101 A010 B010 A101
#  3 A111 B111 A001 B001
#  1 B011 A100 B100 A011
#  0   1         3       5      7
# still gets rotated    \ 4 /
#                               1 \/ 3
#                         /\
#                           / 2 \

Bs           100
55    55    51  55   77   37   77   77    33  77   33   73
010 010 001 010 011 000 011 011 111 011 111 100
7.0,3.0,   3.0,3.0,   7.0,7.0   ,   3.0,3.0,   7.0,7.0,    7.0,7.0,    3.0,7.0,   7.0,7.0,    5.0,5.0   ,    5.0,1.0,   5.0,5.0,   5.0,5.0



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


/*
 * this thing converges really quickly this is more than enough for fp

inline float quicksqrt(float b)
{
    float x = 1.1;
    unsigned char i =0;

    for(;i<16;i++){
        x = (x+(b/x))/2.0;
    }

    return x;
}
*/
/*WARNING: not true euclidean distance
 * compute the distance between two 24 dimensional vectors.
 * The square-root is omitted because the algorithm only needs
 * to know which is closer d(cp, pt.) or d(cp',pt) , for which
 * sqrt(d(cp, pt.)) and sqrt(d(cp', pt.)) inequality holds for positive
 * distances(this is why we keep the squares).
 */
//#define golay

inline float distance(float *cp,float pt[2])
{
 // printf("%f,%f : %f,%f = ",cp[0],cp[1],pt[0],pt[1]);
    float s =(cp[0]-pt[0])*(cp[0]-pt[0]) + (cp[1]-pt[1])*(cp[1]-pt[1]);
   // printf(" %f\n",s);
    return s;

}





/*
 * an integer symbol encoding of an H6 encoder.
 * 0 1 2 3 = 0 1 w w'
 * indexes of the array result in the H6 encoding
 * of a 3 integer symbol character equivalent word.
 * eg: H6CodeWords[0][1][2] = [0,3,2]
 *resulting in the codeword : 0 1 w 0 w' w
 */
static const unsigned char H6CodeWords[4][4][4][3]  = {
    {//0    0       1       w      -w
        {{0,0,0},{1,1,1},{2,2,2},{3,3,3}},//0
        {{1,2,3},{0,3,2},{3,0,1},{2,1,0}},//1
        {{2,3,1},{3,2,0},{0,1,3},{1,0,2}},//w
        {{3,1,2},{2,0,3},{1,3,0},{0,2,1}}},//-w
    {//1    0       1       w      -w
        {{1,3,2},{0,2,3},{3,1,0},{2,0,1}},//0
        {{0,1,1},{1,0,0},{2,3,3},{3,2,2}},//1
        {{3,0,3},{2,1,2},{1,2,1},{0,3,0}},//w
        {{2,2,0},{3,3,1},{0,0,2},{1,1,3}}//-w
    },
    {//w    0       1       w      -w
        {{2,1,3},{3,0,2},{0,3,1},{1,2,0}},//0
        {{3,3,0},{2,2,1},{1,1,2},{0,0,3}},//1
        {{0,2,2},{1,3,3},{2,0,0},{3,1,1}},//w
        {{1,0,1},{0,1,0},{3,2,3},{2,3,2}}//-w
    },
    {//-w   0       1       w      -w
        {{3,2,1},{2,3,0},{1,0,3},{0,1,2}},//0
        {{2,0,2},{3,1,3},{0,2,0},{1,3,1}},//1
        {{1,1,0},{0,0,1},{3,3,2},{2,2,3}},//w
        {{0,3,3},{1,2,2},{2,1,1},{3,0,0}}//-w
    }
};


static const unsigned char H6CodeWordsRev[4][4][4][3]  = {
    {
    //#w    0       1       w      -w
    { {0,0,0}, {3,2,1}, {1,3,2}, {2,1,3} }, //0
    { {2,3,1}, {1,1,0}, {3,0,3}, {0,2,2} },//1
    { {3,1,2}, {0,3,3}, {2,2,0}, {1,0,1} }, //w
    { {1,2,3}, {2,0,2}, {0,1,1}, {3,3,0} }//-w
    },
    {//w    0       1       w      -w
    { {1,1,1}, {2,3,0}, {0,2,3}, {3,0,2} },//0
    { {3,2,0}, {0,0,1}, {2,1,2}, {1,3,3} },//1
    { {2,0,3}, {1,2,2}, {3,3,1}, {0,1,0} },//w
    { {0,3,2}, {3,1,3}, {1,0,0}, {2,2,1} }//-w
    },
    {//w    0       1       w      -w
    { {2,2,2}, {1,0,3}, {3,1,0}, {0,3,1} },//0
    { {0,1,3}, {3,3,2}, {1,2,1}, {2,0,0} },//1
    { {1,3,0}, {2,1,1}, {0,0,2}, {3,2,3} },//w
    { {3,0,1}, {0,2,0}, {2,3,3}, {1,1,2} }//-w
    },
    {//-w   0       1       w      -w
    { {3,3,3}, {0,1,2}, {2,0,1}, {1,2,0} },//0
    { {1,0,2}, {2,2,3}, {0,3,0}, {3,1,1} },//1
    { {0,2,1}, {3,0,0}, {1,1,3}, {2,3,2} },//w
    { {2,1,0}, {1,3,1}, {3,2,2}, {0,0,3} }//-w
    }};



/*
// shaping -.75, -.25,+.25,+.75
//the unit scaled points of 16QAM centered at the origin.
// along with their golay code + parity bit representations
//000, 110 , 001, 111
static const float evenAPts[4][2] = {{-.75, .75},{.25, .75},{.25, -.25},{-.75, -.25}};
//010 100 011 101
static const float oddAPts[4][2]  ={{-.25, .25},{-.25, -.75},{.75, -.75},{.75, .25}};
//000, 110 , 001, 111
static const float evenBPts[4][2] = {{-.25, .75},{.75, .75},{.75, -.25},{-.25, -.25}};
//010 100 011 101
static const float oddBPts[4][2]  = {{.25, .25},{.25, -.75},{-.75, -.75},{-.75, .25}};
*/

//for QAM 64 and whole numbers
//000, 110 , 001, 111
static const float evenAPts[4][2] = {{1.0, 7.0},{5.0, 7.0},{5.0, 3.0},{1.0, 3.0}};
//010 100 011 101
static const float oddAPts[4][2]  ={{3.0, 5.0},{3.0, 1.0},{7.0, 1.0},{7.0, 5.0}};
//000, 110 , 001, 111
static const float evenBPts[4][2] = {{3.0, 7.0},{7.0, 7.0},{7.0, 3.0},{3.0, 3.0}};
//010 100 011 101
static const float oddBPts[4][2]  = {{5.0, 5.0},{5.0, 1.0},{1.0, 1.0},{1.0, 5.0}};


/*
*    this function returns all of the pertinant information
*    from the decoder such as minimum distances, nearest
*    coset leader quadrant, and alternative k-parity distances
*
* these maps are seperated into the quadrants of a cartesian
*  plane now we gotta order these properly
*
* another simple fix is that the quadrants of QAM be abstractly
* defined, and the -,+ of order pairs be used to tile the
* generalized 16bit qam, besides this has to be done anyway
*  so we can get out the real number coordinates in the end
 */
void QAM(float *r, float evenPts[4][2],float oddPts[4][2],float dijs[12][4],float dijks[12][4],unsigned char kparities[12][4]){
//void QAM(float *r, float *evenPts,float *oddPts,float *dijs,float *dijks,int *kparities){


    //the closest even-type Z2 lattice point is used as the
    //coset representatives for all points, not currently used
    //quadrant = [0 for k in range(12)]

    unsigned char i = 0;

    for(;i<12;i++){

        float dist000 = distance(&r[i*2],evenPts[0]);
        float dist110 = distance(&r[i*2],evenPts[1]);
        float dist001 = distance(&r[i*2],evenPts[2]);
        float dist111 = distance(&r[i*2],evenPts[3]);

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
        float dist010 = distance(&r[i*2],oddPts[0]);
        float dist100 = distance(&r[i*2],oddPts[1]);
        float dist011 = distance(&r[i*2],oddPts[2]);
        float dist101 = distance(&r[i*2],oddPts[3]);
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



/*
    computes the Z2 block confidence of the concatonated points projections onto GF4 characters
*/
void blockConf(float dijs[12][4],float muEs[6][4],float muOs[6][4],unsigned char prefRepE[6][4][4],unsigned char prefRepO[6][4][4]){


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


/*here we are looking for the least character in the H6 hexacdoe word
   returns the hexacode word and the wt, for using in locating the least reliable symbol
*/
void constructHexWord(float mus[6][4],unsigned char chars[6],float charwts[6]){

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
    }
}



/*
    this is the minimization over the hexacode function using the 2nd algorithm of  amrani and be'ery ieee may '96
*/
float minH6(unsigned char  y[6],float charwts[6],float mus[6][4]){


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
    unsigned char  candslst[8][6]=  {{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},
                                                        {0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
    //(unsigned char[8][6]) (malloc(8*6*sizeof(unsigned char)));//{
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


        candslst[i+4][0] = H6CodeWordsRev[y[3]][y[4]][y[5]][0];
        candslst[i+4][1] = H6CodeWordsRev[y[3]][y[4]][y[5]][1];
        candslst[i+4][2] = H6CodeWordsRev[y[3]][y[4]][y[5]][2];
        candslst[i+4][3] = y[3] ;
        candslst[i+4][4] = y[4];
        candslst[i+4][5] = y[5] ;



        /*
        candslst[i+4][0] = y[3];
        candslst[i+4][1] = y[4];
        candslst[i+4][2] = y[5];
        candslst[i+4][3] = H6CodeWordsRev[y[3]][y[4]][y[5]][0];
        candslst[i+4][4] = H6CodeWordsRev[y[3]][y[4]][y[5]][1];
        candslst[i+4][5] = H6CodeWordsRev[y[3]][y[4]][y[5]][2];*/
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



/*
    here we are resolving the h-parity. which requires that the overall least significant bit parities equal the
    bit parities of each projected GF4 block. aka column parity must equal 1st row parity
*/
float hparity(float weight,unsigned char hexword[6],unsigned char prefReps[6][4][4],float dijs[12][4],unsigned char oddFlag,unsigned char * codeword){

   unsigned char parity= 0;
   unsigned char i =0;



    for(;i<6;i++){

        //create the golay codeword from the hexacode representation
        codeword[i*4]=prefReps[i][hexword[i]][0];
        codeword[i*4+1]=prefReps[i][hexword[i]][1];
        codeword[i*4+2]=prefReps[i][hexword[i]][2];
        codeword[i*4+3]=prefReps[i][hexword[i]][3];

        //
        parity = parity + prefReps[i][hexword[i]][0];//this should be the highest order bit
    }



    /*this is needed otherwise 2x golay code are emitted*/

    if((parity&1) == oddFlag){
        return weight;
    }


    float leastwt = 1000.0;
    unsigned char least = 0;
    float deltaX;
    unsigned char idx1,idx2;
    i = 0;
    //walk along the codeword again
    for(;i<6;i++){

        idx1 =(codeword[4*i]<<1) +codeword[4*i+1];
        idx2 =(codeword[4*i+2]<<1) +codeword[4*i+3];
        // compute cost of complementing the hexacode representation ^3 of bits
        //select minimal cost complement.
        deltaX = (dijs[2*i][idx1^3] + dijs[2*i+1][idx2^3]) - (dijs[2*i][idx1] + dijs[2*i+1][idx2]);


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

    float least =1000;
    float dif;
    unsigned char argLeast = 0;


    for( ;i <12;i++)
     {
        unsigned char n =(codeword[2*i]<<1)+codeword[2*i+1];
         parity= parity^kparities[i][n];
         codeParity[i] = kparities[i][n];

          dif = dijks[i][n]-dijs[i][n];
          if(dif <= least)
            {
              least = dif;
              argLeast = i;
          }
     }

/*something here as this parity check doesnt fix anything*/
//not sure why this doesnt at least double the set cardinality
    if(parity== Btype ){
        return weight;
    }

    codeParity[argLeast ]=  codeParity[argLeast ] ^1;

    return weight+least;
}





//unsigned char* decode(float r[12][2], float *distance){
unsigned long decodeLeech(float *r, float *distance){

// #####################QAM Dijks ###################
    float* dijs = malloc(sizeof(float)*12*4) ; //{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    float* dijks =malloc(sizeof(float)*12*4) ;// {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    //there is a set for each quarter decoder, and the A/B_ij odd/even
    unsigned char* kparities =malloc(sizeof(unsigned char)*12*4) ;// {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    QAM(r,evenAPts,oddAPts,dijs,dijks,kparities);


    // #####################Block Confidences ###################
    //         0  1    w   W
    float* muEs = malloc(sizeof(float)*6*4*4) ;//{{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}};
    float *muOs = malloc(sizeof(float)*6*4*4) ;//{{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}};
    unsigned char* prefRepE=malloc(sizeof(unsigned char)*6*4*4) ;// = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    unsigned char* prefRepO=malloc(sizeof(unsigned char)*6*4*4) ;// = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};

    blockConf(dijs,muEs,muOs,prefRepE,prefRepO); //just run through both as its faster, but could conserve array allocation

    unsigned char i;

    // #####################Construct Hexacode Word ###################
    unsigned char *y = malloc(sizeof(unsigned char)*6) ;;


    float* charwts = malloc(sizeof(float)*6) ;
    constructHexWord(muEs,y,charwts);


    // #####################Minimize over the Hexacode ###################
    unsigned char* hexword =  malloc(sizeof(unsigned char)*6) ;
    float weight = minH6(y,charwts,muEs);



    //****chars = y = hexword *****
    unsigned char* codeword =  malloc(sizeof(unsigned char)*24);//{0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
    unsigned char* codeParity =  malloc(sizeof(unsigned char)*12) ;

    int winner = 0;

    weight = hparity(weight,y,prefRepE,dijs,0,codeword);//byref

    weight =kparity(weight,codeword,0, codeParity,dijks,dijs,kparities);

    float leastweight = weight;


    //unsigned long leastCodeword;
    unsigned char* leastCodeword = malloc(24*sizeof(unsigned char));
    unsigned long retOpt = 0UL;
    unsigned char b;

    //A is default the least weight decoding
    for(i=0;i<24;i++)
            retOpt = (retOpt)+(codeword[i]<<i);
#ifndef golay
        for(;i<36;i++)
            retOpt = (retOpt)+(codeParity[i-24]<<i);
#endif

    *distance =0;

    //----------------A Odd Quarter Lattice Decoder----------------

    constructHexWord(muOs,y,charwts);;
    weight = minH6(y,charwts,muOs);

    weight = hparity(weight,y,prefRepO,dijs,1,codeword);//byref
    weight =kparity(weight,codeword,0,codeParity,dijks,dijs,kparities);


    if(weight<leastweight)
    {
        leastweight = weight;
        retOpt = 0UL;
        for(i=0;i<24;i++)
                retOpt = (retOpt)+(codeword[i]<<i);
#ifndef golay
            for(;i<36;i++)
                retOpt = (retOpt)+(codeParity[i-24]<<i);
#endif


        *distance =0;
        winner = 1;
    }

    //----------------H_24 Half Lattice Decoder for B points----------------
    QAM(r,evenBPts,oddBPts,dijs,dijks,kparities);
    blockConf(dijs,muEs,muOs,prefRepE,prefRepO);


    //----------------B Even Quarter Lattice Decoder----------------
    constructHexWord(muEs,y,charwts);
    weight = minH6(y,charwts,muEs);

    weight = hparity(weight,y,prefRepE,dijs,0,codeword);//byref
//    printf("BptEvens\n");

    weight =kparity(weight,codeword,1,codeParity,dijks,dijs,kparities);

    if(weight<leastweight){
        retOpt = 0UL;

        leastweight = weight;

        for(i=0;i<24;i++)
                retOpt = (retOpt)+(codeword[i]<<i);
    #ifndef golay
            for(;i<36;i++)
                retOpt = (retOpt)+(codeParity[i-24]<<i);
    #endif

        *distance =1;
        winner = 2;

    }

    //----------------B Odd Quarter Lattice Decoder----------------
    constructHexWord(muOs,y,charwts);
    weight = minH6(y,charwts,muOs);
    weight = hparity(weight,y,prefRepO,dijs,1,codeword);//byref
//    printf("BptOdds\n");

    weight =kparity(weight,codeword,1,codeParity,dijks,dijs,kparities);

    if(weight<leastweight){
        retOpt = 0UL;
        leastweight = weight;
        for(i=0;i<24;i++)
                retOpt = (retOpt)+(codeword[i]<<i);
    #ifndef golay
            for(;i<36;i++)
                retOpt = (retOpt)+(codeParity[i-24]<<i);
    #endif

        *distance =1;
        winner =3;
    }
    *distance = winner;
    //*distance += leastweight;
    free(dijs);
    free(dijks);
    free(kparities);
    free(muEs);
    free(muOs);
    free(prefRepO);
    free(prefRepE);
    free(y);
    free(hexword);
    free(charwts);
    free(codeword);
    free(codeParity);
    free(leastCodeword);

    return retOpt;//leastCodeword;
}
