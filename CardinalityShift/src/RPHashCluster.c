#include "LSHDecoder.c"
#include "leechArrayDecoder.c"
#include "e8Decoder.c"
#include "QAM16Decoder.c"


/*
 * to avoid name conflict, test distance is
 * a euclidean distance function for arbitrary
 * length @len vectors r and q. furthermore
 * unlike@distance() in the main class, it
 * takes the square root.
 */
float testDist(unsigned char * r , float * q,int len ){
  float dist = 0;
  int i;
  for (i=0;i<len;i++)
  {
      //printf("%.1f,",(float)r[i]);
      dist += ((float)r[i]- (float)q[i]) *((float)r[i]- (float)q[i])  ;
  }
  //printf("\n");
return quicksqrt(dist);
}

//more memory conservative 2x time complexity
//main idea is to probe twice, first round maintain counts for all buckets
//second round only store and compute the top k bucket averages
//the complexity is still linear, while memory is O(k*m) instead of O(k*m *||lambda||)
float* rpHash2(float* ret, int numPoints, int dim, Quantizer* q, int hashMod, int cutoff)
{
  int d = q->dimensionality;

    int i,k,j=0;

    int numProbes = 0;
    while((1<<numProbes++) < numPoints  );// poor mans log(numPoints)

    //size of the buckets
    int* buckets = malloc(sizeof(int)*hashMod);



    //sum of decoding distances
    float* distances = malloc(sizeof(float)*hashMod);

    //pass by reference container
    float distance;

    for(;j<numProbes;j++){

      float* M = GenRandomN(d,dim,numPoints);
      for(i=0;i<numPoints;i++)
      {

          long hash = lshHash(&ret[i*dim],dim, 1, hashMod,M, &distance);

          //accumulate hits to this bucket
          buckets[hash]++;
          distances[hash]+=distance;


      }

      free(M);
    }


    //filtering intersections may not be feasible anymore!!!!
    //bucket centroid
    float* bucketsAvgs = malloc(sizeof(float)*k*dim);

    for(j=0;j<numProbes;j++){
      float* M = GenRandomN(d,dim,numPoints);
      for(i=0;i<numPoints;i++)
      {
          long hash = lshHash(&ret[i*dim],dim, 1, hashMod,M, &distance);
    //compute a moving average for current bucket
    //m_{i+1} = (m_{i}+v)/(n+1)

     //if hash is one of the top hashes...
    for(k=0;k<dim;k++)
      bucketsAvgs[hash*dim+k]=(bucketsAvgs[hash*dim+k]*((float)buckets[hash])+ret[i*dim+k])/(((float)buckets[hash])+1.0);
      }
      free(M);
    }



    int *centroidIdx = malloc(sizeof(int)*cutoff);

    //initialize the top cutoff buckets
    for(i=0;i<cutoff;i++)centroidIdx[i] = buckets[i];

    //find the top k
    for(i=cutoff;i<hashMod;i++)
    {
      int argleast = 0;
      int sizeofCandidate = buckets[i];

      for(k=0;k<cutoff;k++)
      {

          //find overlaps in the list k
          if(testDist(&bucketsAvgs[i*dim], &bucketsAvgs[centroidIdx[k]*dim],dim)<2.0){
              //buckets[centroidIdx[argleast]]+= buckets[i];//this needs to be weighted on the two scenarios
              buckets[i]=0;
              k=cutoff;
          }else{

            //checks for biggest buckets
            if(  sizeofCandidate >  buckets[centroidIdx[k]] ){

              //then the kth bucket can be replaced
              //but we need to also keep looking if
              //something is even less

              sizeofCandidate = buckets[centroidIdx[k]];

              //this is the vetor that is getting replaced
              argleast = k;
            }
          }
        }
        //after the least is found it should be in argleast
        if( buckets[i] >buckets[centroidIdx[argleast]] )centroidIdx[argleast] = i;
    }

    float * centroids = malloc(sizeof(float)*dim*k);

    for(k=0;k<cutoff;k++)
    {
      for(j = 0;j<dim;j++)centroids[k*dim+j]=bucketsAvgs[dim*centroidIdx[k] + j];

    }

    free(buckets);
    free(distances);
    free(centroidIdx);
    free(bucketsAvgs);
    return centroids;
}




/*
 * This function projects points until relatively large buckets of points begin to emerge.
 * Serving as an attempt to prove the hypothesis that these points are the approx.
 * density modes of the later discovered clusters.
 * max number of tests, l is the target bucket size before emitting to the parallel system
 * dim is the vector dimensionality
 */
float* rpHash(float* ret, int numPoints, int dim, Quantizer* q, int hashMod, int cutoff)
{
  int d = q->dimensionality;
  int i,k,j=0;

  int numProbes = 0;
  while((1<<numProbes++) < numPoints  );// poor mans log(numPoints)
  //while(numProbes*numProbes<hashMod)numProbes++;<-- Bday Paradox, way to high!


  /*
    <---------------THIS IS THE MAPPER--------------->
   */

  //size of the buckets
  int* buckets = malloc(sizeof(int)*hashMod);
  for(i=0;i<hashMod;i++) buckets[i]=0;
  //bucket centroid
  float* bucketsAvgs = malloc(sizeof(float)*hashMod*dim);
  for(i=0;i<hashMod*dim;i++) bucketsAvgs[i]=0.0f;
  //sum of decoding distances
  float* distances = malloc(sizeof(float)*hashMod);
  //pass by reference container
  float distance;



  long l = time(NULL);
  for(j=0;j<numProbes;j++){
    float* M = GenRandomN(d,dim,numPoints);
    for(i=0;i<numPoints;i++)
    {
        int hash = (int) lshHash(&ret[i*dim],dim, 1, hashMod,M, &distance);

        //accumulate hits to this bucket
        buckets[hash]++;
        distances[hash]+=distance;

        //compute a moving average for current bucket
        //m_{i+1} = (m_{i}+v)/(n+1)
        for(k=0;k<dim;k++)
          bucketsAvgs[hash*dim+k]=(bucketsAvgs[hash*dim+k]*((float)buckets[hash])+ret[i*dim+k])/(((float)buckets[hash])+1.0);
    }
    free(M);
  }

  printf("%u\n",time(NULL)-l);
  int *centroidIdx = malloc(sizeof(int)*cutoff);

  //initialize the top cutoff buckets
  for(i=0;i<cutoff;i++)
      centroidIdx[i] = i;

  //find the top k
  for(i=cutoff;i<hashMod;i++)
  {
    int argleast = 0;
    int sizeofCandidate = buckets[i];
/*
  <---------------THIS IS THE REDUCER--------------->
 */


    for(k=0;k<cutoff;k++)
    {
        //use XOR'd lattice ID as general distance measure
        //find overlaps in the list k
        //can compute as max distance to cluster's average distance                    vvvvv though it may still be gaussian bound apriori, simulate to check.
        if(testDist(&bucketsAvgs[i*dim], &bucketsAvgs[centroidIdx[k]*dim],dim)<quicksqrt(3.0))
          {
            //weighted average nearby buckets
            int s = buckets[centroidIdx[k]] + buckets[i];
#ifdef AVG_CLUSTER
            for(j=0;j<dim;j++)
            {
                bucketsAvgs[centroidIdx[k]*dim+j] = ( (buckets[centroidIdx[k]])*bucketsAvgs[centroidIdx[k]*dim+j]
                                                                                           + (buckets[i]) * bucketsAvgs[i*dim+j]  )/((float)s);
            }
#endif
            buckets[centroidIdx[k]]=s;
            buckets[i] =0;
            k=cutoff;//break
        }else{
          //checks for biggest buckets
          if(  sizeofCandidate >  buckets[centroidIdx[k]] )
            {
            //then the kth bucket can be replaced
            //but we need to also keep looking if
            //something is even less
            sizeofCandidate = buckets[centroidIdx[k]];

            //this is the vetor that is getting replaced
            argleast = k;
          }
        }
      }

      //after the least is found it should be in argleast
      if( buckets[i] >buckets[centroidIdx[argleast]] )centroidIdx[argleast] = i;
  }

  //load return vector of centroids
  float * centroids = malloc(sizeof(float)*dim*k);

  for(k=0;k<cutoff;k++)
  {
    for(j = 0;j<dim;j++)centroids[k*dim+j]=bucketsAvgs[dim*centroidIdx[k] + j];

    printf("%i:%i\n",centroidIdx[k],buckets[centroidIdx[k]]);
  }

  free(buckets);
  free(distances);
  free(centroidIdx);
  free(bucketsAvgs);
  return centroids;

}
