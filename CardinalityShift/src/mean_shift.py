
H  = 1; #if useKNNToGetH is false, H will be used else knn will be used to determine H
threshold = 1e-3;
numPoints = 300;
k = 100; % when we use knn this is the k that is used.

numClasses = 0;
pointsToCluster2 = [];
centroids = [];

def getDataPoints(numPoints,numDimensions):
        numClasses = randi([2,8],1) # generate a random number between 2 and 8
        dataPoints = []
        curNumberOfPoints = 0

        randomMeans = []
        for i in range(numClasses):
                randomMean = randi([0,100], 1) * rand(1)
                randomStd = randi([1,2]) * rand(1)
                curGaussianModelPoints = randomMean + randomStd .* randn(numPoints,numDimensions)
                dataPoints = [dataPoints;curGaussianModelPoints]
                randomMeans = [randomMeans,randomMean]
                pointsToCluster2(curNumberOfPoints +1 : (curNumberOfPoints+numPoints)) = i
                curNumberOfPoints = curNumberOfPoints + numPoints
                centroids = randomMeans
                
        return dataPoints

def sqdist(a,b):
        #taken from demo code of class
        aa = sum(a.*a,1)
        bb = sum(b.*b,1)
        ab = a.T*b
        return abs(repmat(aa.T,[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab)

def getBandWith(curDataPoint,origDataPoints,euclideanDist,useKNNToGetH):
    if (useKNNToGetH == false):
        bandwidth = H
    else:
        [sortedEuclideanDist,indexInOrig] = sort(euclideanDist)
        bandwidth = sortedEuclideanDist(k)
    return bandwidth

def getClusters(origDataPoints,finalDataPoints):
    [numSamples,numFeatures] = len(finalDataPoints),len(finalDataPoints[0])
    clusterCentroids  = []
    pointsToCluster = []
    numClusters = 0

    for i in range(numSamples):
        closestCentroid = 0 
        curDataPoint = finalDataPoints[i]

        for j in range(numClusters):
            distToCentroid = sqrt(sqdist(curDataPoint.T,clusterCentroids[j].T))
            if (distToCentroid <  4 * H):
                closestCentroid = j
                break

        if (closestCentroid > 0):
            pointsToCluster[i] = closestCentroid
            clusterCentroids[closestCentroid] = 0.5 * (curDataPoint + clusterCentroids[closestCentroid])
        else:
            numClusters = numClusters + 1
            clusterCentroids[numClusters] = finalDataPoints[i]
            pointsToCluster[i] = numClusters
 
    return clusterCentroids,pointsToCluster


def plotPoints(clusterCentroids,pointsToCluster,origDataPoints):
    from pylab import plot
    numSamples,numFeatures = len(origDataPoints),len(origDataPoints[0])

    if (numFeatures == 2):
        # plot the original Data Points
        h = figure(1);
        #hold on;
        for i in range(numClasses):
            colourIndex = mod(i,sizeOfColors) + 1
            allElemsInThisCluster = find(pointsToCluster2 == i)
            allOrigPointsInThisCluster = origDataPoints[allElemsInThisCluster]
            plot([ i[1] for i in allOrigPointsInThisCluster ],[i[2] for i in allOrigPointsInThisCluster])
    
        plot(centroids,centroids)


        h = figure(2)
        #plot the original Data Points
        for i in range(size(clusterCentroids,1)):
            colourIndex = mod(i,sizeOfColors) + 1 
            allElemsInThisCluster = filter(i,pointsToCluster)
            allOrigPointsInThisCluster = origDataPoints[allElemsInThisCluster]
            plot([ i[1] for i in allOrigPointsInThisCluster ],[i[2] for i in allOrigPointsInThisCluster]



from math import exp
def doMeanShift(dataPoints,useKNNToGetH):
    numSamples,numFeatures = len[dataPoints],len[dataPoints[0]]
    origDataPoints = dataPoints
    for i in range(numSamples):
        diffBetweenIterations = 10
        while (diffBetweenIterations > threshold):
                curDataPoint = dataPoints[i]
                euclideanDist = sqdist(curDataPoint.T,origDataPoints.T)
                bandwidth = getBandWith(origDataPoints[i],origDataPoints,euclideanDist,useKNNToGetH)
                kernelDist = [exp(-euclideanDist[i]/bandwidth[i]^2) for i in xrange(len(euclideanDist))]
                numerator = [kernelDist[i] * origDataPoints[i] for i in xrange(len(kernelDist))]  
                denominator = sum(kernelDist)
                newDataPoint = [numerator[i]/denominator[i] for i in xrange(len(numerator))]
                dataPoints[i] = newDataPoint
                diffBetweenIterations = abs(curDataPoint - newDataPoint) 

    clusterCentroids,pointsToClusters = getClusters(origDataPoints,dataPoints)
    plotPoints(clusterCentroids,pointsToClusters,origDataPoints)
    
    return origDataPoints,dataPoints
