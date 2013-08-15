Copyright 2013 Lee Carraher. All rights reserved.

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



_Intuition_:
------------
The basic principle is random projection and discrete space quantization.
Multiple Gaussian random projections allow every vectors to occupy a fuzzy
region around their actual location. The discrete space quantizers hard margin
boundaries will only decode a point if it is specifically in its error
radius. But by
allowing the point to occupy a fuzzy region around its actual location, we get
a probability of being in- or out-of-a-bucket that correlates with the vectors'
real distance from the center of the space quantizer. By doing this a bunch of
times, we can get a sort of smear of all the vectors. And buckets that
on average
generate more hits, tend to be ones near large populations of vectors. And
ff there are large populations of vectors, we can consider the average of
those vectors in that population to be a pretty good centroids. The
neat thing is-is
that we don't have to communicate vectors between nodes. we just have to
communicate the bucket counts.

    
    
_prerequisites_:
----------------
1. n vectors of d-dimensional data
2. cluster of compute nodes each containing part of the n vectors
(does not necessarily need to be homogenized)
3. standard map reduce or mpi network routing (just need gather and distribute)

_algorithm_:
------------
1. Each node generates log(n) random gaussian projection matrices, M_d,m
* where m corresponds to the quantizer subspace (leech - 24, e8 - 8, BW - 32...)
* we choose log(n) times as this flattens the gaussians probability
density function.
2.  Each node projects its vectors log(n)-times using the log(n) M_d,m  matrices
3. Each node computes the space quantizer decodings for every
projection of the vectors it has
4. Store the decoding with the vectors
5. Count the occurrences of all decodings
6. Gather the {decoding, count}-pairs for all nodes
7. Choose the top eps-k larges {decoding,count}-pairs
8. Distribute the k larges decodings to every node
9. Each node applies the decodings as classification sets to the vectors it has
10. Each node computes an average centroid for the classifications sets
11. The average centroids are gathered and averaged over all nodes
12. Within the eps-k centroids, choose the top k that do not overlap
*this last step needs work but my not be needed
END

_analysis_
------------
Matrix generation: `O(log(n)* d*m)`  
Projection: `O(log(n) * n*d*m)`  
Decoding: `O(c*log(n)*n)`  
Count: `O(log(n)*n)`  
Gather: `O(1)`  
Distribute: `O(1)`  
Average: `O(log(n)*n)`  
Gather: `O(1)`  
*projection dominates complexity so overall  
      `O(log(n) * n*d*m)`  
* and this parallelizes over the number of nodes well  


Test by compiling and running LSHTests.c
----------------------------------------
USAGE:
------
1. Run a random example with 20 cluster of 3000 x 1000 dimensional vectors
>  ./rptest 

2. Generate random data file with 20 cluster of numberOfPtsPerCluster x 1000 
dimensional vectors
>  ./rptest numberOfPtsPerCluster
outputs two files:  
ORIGDATA.mat - the original cluster centers  
TEMPOUT.mat  - all data vectors  

3. Find K clusters in a *.mat(see tinyBLAS for format) datafile. Compare with
ORIGDATA.mat if available, and use a HASHMOD buckets.
>  ./rptest INPUTFILE.mat K [ORIGDATA.mat] [HASHMOD]
 
OUTPUT:
-------

`Leech Dec : First 5 pts`  
3672517440 : 0.29,0.61,-0.31,0.06,0.93  
732610504 : -0.35,0.25,0.23,-0.24,0.35  
2701931211 : 0.98,-0.66,0.51,0.83,-0.58  
3015201805 : 0.82,0.13,0.07,0.80,0.81  
963235791 : -0.94,0.65,-0.46,0.62,-0.61  
1179984320 : 0.42,0.15,0.58,-0.35,-0.77  
3067812570 : -1.00,0.38,-0.23,0.27,-0.06  
2321353263 : 0.91,0.92,0.37,0.56,0.59  
174042856 : -0.83,-0.64,-0.45,0.75,0.84  
2306016502 : -0.72,-0.41,0.22,0.55,0.67  
1288249006 : 0.54,0.35,0.86,-0.49,0.87  
1036683323 : 0.86,0.21,0.38,0.97,-0.53  
1252308019 : 0.34,-0.89,0.36,0.10,0.81  
924406298 : -0.54,-0.01,-0.05,-0.94,0.64  
4142256666 : 0.07,-0.90,0.22,-0.56,0.60  
364237980 : -0.39,0.15,0.20,-0.38,0.51  
1116141030 : -0.05,0.73,-0.58,-0.11,0.03  
3801150847 : -0.16,-0.12,-0.17,-0.84,0.60  
1388838905 : -0.55,0.09,0.44,0.24,-0.37  
952430016 : 0.24,0.99,-0.84,0.92,-0.92  
  
180   `<------- how long it took>`  
3607:2011  
4243:2150  
2676:4282  
138:4076  
432:5849  
2078:2844  
5187:1340  
1311:1929      `<---- Bucket labels and sizes>`  
6820:1387  
297:5091  
105:8368  
2139:3529  
46:4091  
262:3361  
7901:1439  
935:5105  
1209:5235  
479:3183  
5909:3953  
2931:2063  
  
`Leech Dec : NN : Dist to NN : First 5 pts`  
1388838905 : 18 : 1.2389 : -0.53,0.01,0.40,0.23,-0.38  
1749329861 : 19 : 4.5600 : 0.07,0.91,-0.78,0.83,-0.76  
2309948616 : 9 : 2.3422 : -0.69,-0.34,0.21,0.46,0.68  
1978503838 : 14 : 1.0998 : -0.09,-0.87,0.22,-0.48,0.57  
1288249006 : 10 : 1.6497 : 0.45,0.37,0.92,-0.56,0.84  
3155182681 : 2 : 2.4270 : 0.79,-0.66,0.38,0.59,-0.31  
1973894816 : 1 : 4.1429 : -0.27,0.38,0.07,-0.01,0.15  
1059323903 : 4 : 9.0860 : -0.64,0.38,-0.20,0.37,-0.31  
952430016 : 19 : 0.9352 : 0.25,0.97,-0.84,0.90,-0.89  
1037493771 : 11 : 2.6351 : 0.72,0.16,0.30,0.80,-0.40  
2827277890 : 0 : 2.4022 : 0.40,0.51,-0.11,0.11,0.79  
174042856 : 8 : 1.1692 : -0.79,-0.63,-0.41,0.69,0.89  
289621565 : 1 : 1.1929 : -0.36,0.20,0.26,-0.44,0.36  
174042856 : 8 : 2.7615 : -0.85,-0.69,-0.35,0.63,0.58  
3141811937 : 2 : 0.3774 : 0.97,-0.65,0.50,0.86,-0.56  
1745668655 : 12 : 1.1327 : 0.37,-0.89,0.28,0.06,0.80  
3052918175 : 15 : 2.2001 : -0.44,0.09,0.13,-0.25,0.61  
1340400192 : 0 : 0.8403 : 0.26,0.70,-0.22,0.04,0.93  
3067812570 : 6 : 0.5696 : -1.00,0.36,-0.22,0.27,0.04  
3801150847 : 17 : 5.4913 : -0.19,-0.04,-0.10,-0.57,0.41  
  
Results:
--------
clustered 1 million 30 dimensional vectors, in under 29 minutes, with average deviation from true gaussian centers of: 0.01089
and average cluster distance : 7.416
  
appears to be linear in `O(n)`
  
> numVectors    = [1000,10000,20000,50000,100000,500000,1000000]
> timeToCluster = [1.888,17.709,36.678,88.226,173.415,885.527,1727.209]
> ratios =>  [0.001888, 0.001771, 0.001834, 0.001765, 0.001734, 0.001771, 0.001727]
  

> BucketSizes = 1015,1018,1018,1053,1069,1080,1133,1133,1145,1147,1166,1169,1176,1243,1252,1265,1279,1282,1284,1314,1340,1387,1396,1439,1456,1463,1465,1493,1535,1590,1637,1655,1683,1708,1710,1748,1756,1768,1782,1820,1846,1859,1862,1871,1909,1929,2011,2012,2022,2035,2063,2068,2085,2090,2090,2096,2120,2140,2150,2169,2253,2267,2309,2552,2566,2655,2781,2796,2844,2863,2865,3001,3085,3128,3172,3183,3202,3241,3326,3355,3361,3435,3461,3473,3506,3529,3562,3565,3642,3653,3706,3713,3741,3805,3818,3953,3954,3965,3997,4032,4076,4091,4199,4282,4393,4430,4437,4464,4516,4518,4563,4696,4723,4742,4867,4879,5091,5105,5199,5221,5235,5253,5338,5384,5403,5414,5439,5441,5442,5551,5593,5672,5676,5743,5807,5849,5958,6051,6067,6111,6159,6449,6618,6881,6911,7030,7284,7539,7549,7631,7686,7853,8296,8368,9081,9085,9411,9648,9864,11693
> DistanceToRealCenters =
6.7981,2.9852,11.1027,2.9855,7.9459,4.9635,6.141,4.2617,5.0781,10.6027,12.6425,7.1342,12.693,3.2326,10.5724,3.99,5.9386,12.7211,0.5472,7.216,4.1429,0.9352,3.3193,0.3774,8.0976,0.9286,2.8413,11.5489,3.6384,7.856,11.7422,3.1929,10.0629,3.7097,4.5633,11.1915,5.9227,5.616,2.9676,1.7073,0.999,13.5599,1.188,3.3922,9.7363,9.086,1.2389,2.8102,4.5541,2.056,5.4913,0.7376,12.0435,1.2015,3.4305,3.3493,1.8083,7.9469,4.56,2.2245,0.9538,8.1558,3.3112,2.4028,1.4665,2.8312,2.5354,2.8253,2.427,0.7663,0.9066,0.4932,0.9127,1.9151,0.5948,0.8403,1.0999,2.1139,1.166,0.9306,2.7615,2.5642,3.3439,3.0404,0.9467,1.1692,0.6046,0.6552,1.564,1.0078,1.3472,0.682,1.9528,0.7839,2.3456,0.5696,0.6821,0.6665,0.4977,0.258,1.0998,1.1929,1.0726,2.3422,2.4371,0.4117,1.4277,1.6005,1.2011,0.5797,0.6777,1.2681,0.8703,2.115,0.9287,0.6941,2.6351,1.1327,1.4368,1.7301,2.2001,1.1894,2.1747,1.4104,2.2489,0.8898,0.7987,2.5369,1.2419,2.3636,1.1295,1.992,0.8772,1.7901,2.12,1.6497,1.2524,0.5821,0.7416,0.6894,1.5614,1.443,1.7759,2.5949,2.4267,2.0128,1.5872,2.5758,0.8316,1.5234,1.6482,1.145,2.1046,2.4022,1.0958,0.9588,2.074,1.4716,1.22,0.5489



