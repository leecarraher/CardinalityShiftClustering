
/*
 * Generate a 'good enough' gaussian random variate.
 * based on central limit thm , this is used if better than
 * achipolis projection is needed
 */
double randn() {
  int i;
  float s = 0.0;
  for(i = 0;i<6; i++)s+=((float)rand())/RAND_MAX;
  return s - 3.0;

}


/*
 *print @size values of an integer vector @v
 */
inline void printVecI(unsigned char* v,int size){
  unsigned char i = 0;
  for(;i<size;i++){
      printf("%i ",v[i]);
      if(i<size-1)printf(",",v[i]);
  }

  printf("\n");
}

/*
 *print @size values of an float vector @v
 */
inline void printVecF(float* v,int size){
unsigned char i = 0;
for(;i<size;i++){
    printf("%.2f",v[i]);
    if(i<size-1)printf(",",v[i]);
}
printf("\n");
}


/*
 * print the binary expansion of a long integer. output ct
 * sets of grsize bits. For golay and hexacode set grsize
 * to 4 bits, for full leech decoding with parity, set grsize
 * to 3 bits
 */
static void print(unsigned long ret,int ct,int grsize){
    int i,j;
    for(i=0;i<ct;i++)
    {
        for(j=0;j<grsize;j++)
        {
            printf("%d",ret&1);
            ret=ret>>1;
        }
        printf(" ");
    } printf("\n");
}

int weight(unsigned long u){
  int ret = 0;
  while(u>0)
   {
      if(u&1 == 1)ret++;
      u = (u>>1);
  }
  return ret;
}





void convertToCoords(unsigned long c,unsigned char Bpoint, unsigned char* point)
{
 /*
#  7 A000 B000 A110 B110
#  5 B101 A010 B010 A101
#  3 A111 B111 A001 B001
#  1 B011 A100 B100 A011
#  0   1         3       5      7
*/
  // since a -1 to 1 scaled 16 QAM is used the below
  // integer scaling is not needed, but is useful for
  // visualizing so it remains in test.
  //    01    23   45   67     89  10 11      12 13  14 15
  //    000  001 010  011 100 101         110      111
  unsigned char axCoords[] = {1,5, 3,7,3,7,5,1 };
  unsigned char ayCoords[] = {7,3,5,1,1,5,7,3};
  unsigned char bxCoords[] = {3,7,5,1,5,1,7,3};
  unsigned char byCoords[] = {7,3, 5,1,1,5,7,3};


  int i = 11;


  if(Bpoint==0){
      for(;i>-1;i--){
            point[i*2]= axCoords[c&7];
            point[i*2+1]=ayCoords[c&7];
            c = c>>3;
      }
  }

  for(;i>-1;i--){
        point[i*2]= bxCoords[c&7];
        point[i*2+1]=byCoords[c&7];
        c = c>>3;
  }

}


/*
 * Create a random float vector @r of length @len and
 * @sparseness
 */
void genRandomVector(int len,float sparesness,float* r)
{
  while(len>0){
      if(((float)rand())/((float)RAND_MAX)< sparesness)
          r[len-1]=2.0*((float)rand())/RAND_MAX-1.0;
      else r[len-1]=0.0;
      len--;
  }
}




/*
 * Lets prove we are not accidentally finding clusters.
 */
void shuffle(float* M,int dim,int len)
{
  int j,i = 0;
  for(;i<len;i++)
  {
      int r = rand()/RAND_MAX;
      for(j = 0;j<dim;j++)
      {
          float temp = M[r*dim+j];
          M[r*dim+j] = M[i*dim+j];
          M[i*dim+j] = temp;

         // xor swap
        /*
          M[r*dim+j] ^= M[i*dim+j] ;
          M[i*dim+j] ^= M[r*dim+j];
          M[r*dim+j] ^= M[i*dim+j] ;
          */
      }

  }

}

/*
 * Create a scaled histogram for the counts data
 * inx and the indexed data iny. b is the number of
 * buckets, and len is the length of iny and inx.
 */

float histogram(int b, int len, float* inx,float* iny)
{
  int i,j;
  float max = -RAND_MAX;
  float min = RAND_MAX;//big number

  int* cts = malloc(sizeof(int)*b);
  float * out =malloc(sizeof(float)*b);

  //initialize output
  for (i=0;i<b; i++)
   {
      cts[i]=0;
      out[i]=0;
   }


  for (i=0;i<len; i++)
   {
      if(inx[i]>max)max =(float) inx[i];
      if(inx[i]<min)min =(float) inx[i];
   }
  float wsize = ((float)max-min) /((float) b);

  for (i=0;i<len; i++)
  {
      j = 0;
      while(j+1<b && inx[i] > min+wsize*j)j++;
      out[j] += (float)iny[i];
      cts[j]++;
  }

  //scale averages
  for (i=0;i<b; i++)
      if(cts[i]!=0)out[i]=out[i]/((float)cts[i]);
  printVecF(out,b);
  for (i=0;i<b; i++)printf("%.0f, ",wsize*i);printf("\n");
  return wsize;
}


int listSearchAndAdd(unsigned long query, unsigned long* list,int len)
{
  list[len]=query;
  while(--len)
      if(list[len]==query)
        return len;

  return len;
}



/*
 * Return the index of the nearest neighbor to vector v
 */
int NN(float* v, float*M,int dim,int len)
{

  int i,j = 0;
  float dist = 0.0f;
  int argmin = 0;
  float min = 10000.0f;
  float* q;
  for(;j<len;j++)
    {
      dist =0.0f;
      q = &M[j*dim];
      for (i=0;i<dim;i++)
      {
          dist += (v[i]- q[i]) *(v[i]- q[i])  ;
      }
      if(dist<min){
          min = dist;
          argmin = j;
      }

  }
  return argmin;
}


/*
 * generate random centers
 */
float * generateRandomCenters(int d,int clu)
{
  int i;
  float *clusterCenters =(float*) malloc(sizeof(float)*d*clu);
  for (i=0;i<clu;i++){
      genRandomVector(d,1.0,&clusterCenters[i*d]);
  }

  return clusterCenters;
}



/*
 * Very optimized Gaussian random cluster generator
 */
float* generateGaussianClusters(int part,int d,int clu,float* clusterCenters){//, float *clusterCenters){//, float *ret){

    float *ret = (float*)malloc(sizeof(float*)*d*part*clu);

    int i,j,k;

    for (i=0;i<clu;i++){
        //must be positive and near 1, around .1 , higher numbers results in data out of range ~(-1:1)
        float variance = .1;//*((float)rand())/RAND_MAX;
        float* ct = &clusterCenters[i*d];
        //these are the partitions
        for(j = 0;j<part;j++)
          {
            //these are the individual vectors
            for( k=0;k<d;k++)
              {
                float pts = (randn()*4.0);
                ret[i*part*d +j*d+k]= (pts*variance+ ct[k]);
            }

        }
    }

    return ret;
}

