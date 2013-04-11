from pylab import *
from numpy import *

def GMM():

    """ Fits two Gaussians to data using the EM algorithm """
    N = 1000
    #ion()
    
    y = 1.*zeros(N)
    # Set up data
    out1 = random.normal(6,1,N)
    out2 = random.normal(3,1,N)
    out3 = random.normal(3,5,N)
    choice = random.rand(N)
    
    w = [choice>=0.5]
    y[w] = out1[w]	
    w = [choice<0.5]
    y[w] = out2[w]
    w = [choice<0.5]
    y[w] = out3[w]
    
    clf()
    hist(y,fc='0.5')
    
	# Now do some learning

	# Initialisation
    mu1 = y[random.randint(0,N-1,1)]
    mu2 = y[random.randint(0,N-1,1)]
    mu3 = y[random.randint(0,N-1,1)]
    s1 = sum((y-mean(y))**2)/N
    s2 = s1
    s3 = s1
    pi = 0.5

	# EM loop
    count = 0
    gamma = 1.*zeros(N)
    nits = 300

    ll = 1.*zeros(nits)
	
    while count<nits:
        count = count + 1

    	# E-step
        for i in range(N):
            gamma[i] = pi*exp(-(y[i]-mu2)**2/s2) / ((1-pi) * exp(-(y[i]-mu1)**2/s1) +  pi* exp(-(y[i]-mu3)**2/s3))
        
    	# M-step
        mu1 = sum((1-gamma)*y)/sum(1-gamma)
        mu2 = sum(gamma*y)/sum(gamma)
        mu3 = sum(gamma*y)/sum(gamma)
        s1 = sum((1-gamma)*(y-mu1)**2)/sum(1-gamma)
        s2 = sum(gamma*(y-mu2)**2)/sum(gamma)
        s3 = sum(gamma*(y-mu3)**2)/sum(gamma)
        pi = sum(gamma)/N
        	
        
    x = arange(-2,2,0.01)
    y = 350*pi*exp(-(x-mu1)**2/s1) + 350*(1-pi)*exp(-(x-mu2)**2/s2) + 350*(1-pi)*exp(-(x-mu3)**2/s3)
    plot(x,y,'k',linewidth=2)
    figure(), plot(ll,'ko-')
    show()
    
GMM()
