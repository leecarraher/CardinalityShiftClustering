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
The basic principle is random projection and discrete space quantization.
Multiple gaussian random projections allow every vectors to occupy a fuzzy
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
that we dont have to communicate vectors between nodes. we just have to
communicate the bucket counts.


_prerequisites_:
1. n vectors of d-dimensional data
2. cluster of compute nodes each containing part of the n vectors
(does not necessarily need to be homogenized)
3. standard map reduce or mpi network routing (just need gather and distribute)

_algorithm_:
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
M generation: O(log(n)* d*m)
projection: O(log(n) * n*d*m)
decoding: O(c*log(n)*n)
count: O(log(n)*n)
gather: O(1) * actually log(number of nodes, but comparatively small)
distribute: O(1) * "
             "
average: O(log(n)*n)
gather:  O(1)


projection dominates complexity
so overall
 O(log(n) * n*d*m) , and this parallelizes over the number of nodes well



Test by compiling and running LSHTests.c

it will generate 5 clusters:

-0.26,-0.73,0.17,0.47,0.36,-0.95,0.45,-0.86,-0.90,-0.79,0.59,0.08,-0.01,-0.64,0.09,-0.33,-0.39,-0.97,-0.57,-0.96,-0.28,0.33,0.02,-0.20
-0.44,0.16,0.57,-0.06,-0.43,0.54,-0.04,-0.65,0.73,-0.44,0.75,-0.47,0.60,-0.23,0.35,0.46,0.12,0.73,-0.76,0.90,0.18,0.56,-0.59,-0.98
-0.59,-0.62,0.56,-0.12,-0.31,0.94,-0.50,-0.14,-0.95,0.13,-0.61,-0.75,-0.15,-0.23,-0.14,0.84,0.65,-0.42,0.97,-0.47,-0.69,-0.35,-0.85,0.29
0.56,-0.69,0.16,0.35,0.50,-0.82,-0.69,0.84,-0.96,0.83,-0.83,0.42,0.05,0.80,-0.95,-0.38,-0.50,0.66,-0.62,-0.08,-0.67,-0.30,0.63,0.84
-0.22,-0.38,0.24,0.09,0.38,0.03,-0.42,0.23,-0.05,0.43,-0.95,0.70,-0.40,0.66,-0.69,0.31,-0.29,0.91,0.59,0.18,0.17,-0.64,0.51,-0.3

then output what CSC thinks the clusters should be:


vv----------these three seemed like they can be merged---------vv
0.47,-0.63,0.20,0.21,0.46,-0.77,-0.79,0.87,-0.96,0.89,-0.84,0.48,0.15,0.80,-0.94,-0.38,-0.49,0.72,-0.45,-0.13,-0.76,-0.32,0.53,0.86
0.66,-0.74,-0.03,0.35,0.57,-0.87,-0.67,0.74,-0.91,0.88,-0.91,0.50,0.09,0.83,-0.99,-0.37,-0.48,0.55,-0.57,-0.01,-0.61,-0.39,0.62,0.78
0.66,-0.58,0.20,0.32,0.47,-0.90,-0.64,0.86,-0.98,0.80,-0.78,0.50,-0.02,0.83,-0.97,-0.47,-0.48,0.69,-0.67,-0.05,-0.73,-0.24,0.66,0.83


vv----------these are all the same----------vv but average out to the
one that starts with -0.22
-0.22,-0.38,0.24,0.10,0.38,0.03,-0.42,0.23,-0.05,0.43,-0.95,0.70,-0.40,0.65,-0.69,0.32,-0.29,0.91,0.59,0.18,0.17,-0.64,0.51,-0.34
-0.22,-0.38,0.24,0.09,0.38,0.03,-0.42,0.23,-0.05,0.43,-0.95,0.70,-0.40,0.66,-0.69,0.31,-0.29,0.91,0.59,0.18,0.17,-0.64,0.51,-0.34
-0.22,-0.38,0.24,0.09,0.38,0.03,-0.42,0.23,-0.05,0.43,-0.95,0.70,-0.40,0.66,-0.69,0.31,-0.29,0.91,0.59,0.18,0.17,-0.64,0.50,-0.33
^^----------these are all the same----------^^

vv ---------- these are pretty good ----------vv
-0.59,-0.62,0.55,-0.12,-0.31,0.94,-0.51,-0.14,-0.95,0.13,-0.61,-0.75,-0.15,-0.23,-0.14,0.85,0.65,-0.42,0.97,-0.47,-0.69,-0.35,-0.85,0.29
-0.44,0.18,0.49,-0.12,-0.31,0.56,-0.03,-0.72,0.65,-0.42,0.77,-0.54,0.67,-0.11,0.30,0.55,0.06,0.70,-0.77,0.95,0.18,0.59,-0.53,-0.98

vv----------ugly duckling----------vv
-0.55,0.07,0.68,0.01,-0.52,0.65,-0.11,-0.72,0.65,-0.42,0.82,-0.52,0.56,-0.22,0.40,0.43,0.14,0.78,-0.73,0.80,0.21,0.63,-0.59,-1.03


