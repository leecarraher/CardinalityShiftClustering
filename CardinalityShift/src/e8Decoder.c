
/*WARNING: not true euclidean distance
 * compute the distance between two 24 dimensional vectors.
 * The square-root is omitted because the algorithm only needs
 * to know which is closer d(cp, pt.) or d(cp',pt) , for which
 * sqrt(d(cp, pt.)) and sqrt(d(cp', pt.)) inequality holds for positive
 * distances(this is why we keep the squares).
 */
inline float distance8(float cp[8],float pt[8])
{
    float s =(cp[0]-pt[0])*(cp[0]-pt[0])
    + (cp[1]-pt[1])*(cp[1]-pt[1])
    + (cp[2]-pt[2])*(cp[2]-pt[2])
    + (cp[3]-pt[3])*(cp[3]-pt[3])
    + (cp[4]-pt[4])*(cp[4]-pt[4])
    + (cp[5]-pt[5])*(cp[5]-pt[5])
    + (cp[6]-pt[6])*(cp[6]-pt[6])
    + (cp[7]-pt[7])*(cp[7]-pt[7]);

    return s;

}


void decodeD8(float* r,float* rInt){
    float dist = 0.0;
    float * rDist = malloc(sizeof(float)*8);;

    unsigned char i;
    for( i=0; i<8;i++)
      {
        unsigned char tr = (unsigned char)r[i];
        if (r[i]>tr+.5){
            rInt[i]=(float)tr+1.F;
            rDist[i] = (tr+1)-r[i];
        }
        else{
            rInt[i]=(float)tr;
            rDist[i] = r[i]-tr;
        }

    }


    unsigned char sum = 0;
    unsigned char argmax = 0;
    float max = 0;
    for( i=0; i<8;i++)
      {
        if(rDist[i]>max)
          {
            max = rDist[i];
            argmax = i;
        }
        sum+=rInt[i];
    }

    if (sum&1)
      {
        if (rInt[argmax]>r[argmax]) rInt[argmax]=rInt[argmax]-1;
        else rInt[argmax]=rInt[argmax]+1;
    }
    free(rDist);

}

unsigned int decodeE8(float* r,float* dist){
    /*
        this decoder uses the D8 partition to decode E8 lattice
        it returns an integer label for the lattice point
*/


    float* nr = malloc(8*sizeof(float));
    float * hr = malloc(8*sizeof(float));
    float * rshift = malloc(8*sizeof(float));;
    unsigned char i;

    //very simple scaling from [-1,1] to [0,2]
    for( i=0;i<8;i++)
        r[i] =r[i]+1.0;

    decodeD8(r,nr);

    for( i=0;i<8;i++)
        rshift[i] =((float)r[i] -.5);

    decodeD8(rshift,hr);


    for( i=0;i<8;i++)rshift[i] =(float)hr[i] +.5;
    free(hr);

     float nrd = distance8(nr,r);
     float hrd =distance8(rshift,r);



     if (hrd<nrd){

         free(nr);
         nr = rshift;
         *dist += hrd;
     }
     else {
         *dist += nrd;
         free(rshift);
     }

     long s = 0;
    //create an int representation
    //0->00 : .5->01 :1.->10 : 1.5->11 2.->100 2.5
    for( i=0;i<8;i++){
        s<<=3;
        s+=((int)(nr[i]*2));
    }

    free(nr);


    return s;

}


