'''
Created on Apr 24, 2013

@author: lee
'''

import genGauss
import crKNN


if __name__ == '__main__':
    pass

    data, cluCntrs = genGauss.getDataPoints(100, 24, 10)
    DBs,Ms = crKNN.buildDB(data,-1.0,11.0)

    q = [1,2,4,2,1,4,2,1,3,4,2,1,4,3,2,1,2,1,1,4,5,5,1,1]
    print len(crKNN.query(q,DBs,Ms,len(data),10,-1.0,8.0))