def mult(n):
    guess = int(n**.5)
    while not (n%guess ==0):
        guess = guess +1
    return n/guess,guess
    
def nextPow2(n):
    if not (n & (n - 1)) == 0:
        n=n-1
        n |= n >> 1
        n |= n >> 2
        n |= n >> 4
        n |= n >> 8
        n |= n >> 16
        n=n+1
    return n

def commGen(ID,acc_mod,pivot):
    '''
        Generative method for calculating non-overlapping communication
        partners per round in a 2^n-node network.
        
        This is the function each node would call per round.
        due to the equal spliting of a bitstring at each index
        the only synchronization that need occur is at the
        step level, all other nodes are in listening mode
    '''
    pivot_mask = 1<<(pivot-1)
    if (ID&pivot_mask) > 0:
        cp = ID & (pivot_mask-1)
        acc = (ID>>pivot)
        t = []
        i = 0

        while i<acc_mod:
            #this is where we would communicate with the append node
            t.append(   ( ((acc+i)&(acc_mod-1))<<pivot )+cp)
            '''
            # here's whats going on:
            # first increment acc by 1
            # n& (acc_mod-1) is equivalent to n % acc_mod
            # shift 1 bit and put the complement of the current
            # bit there '0' which happens intrinsically with a shift
            # now shift some space to append the lower nibble to the bitstring)
            '''
            
            i=i+1
        return t 

 
def generalRun(n,fnc,vargs):
    '''
        this is the sequential run version of the comm function.
        in general the while loop would synchronize nodes and the
        for ID loop would be done in parallel on the |nodes|
    '''
    mod=n>>1   
    pivot = 1
    '''
    # start the pivot at the first bit,
    #mod says there will be n/2 transfers for the first round
    #(n/2)/2 for the second and (n/2)/2 .../2 subsequently
    '''
    while mod>0:
        actives = []
        passives = []
        for ID in xrange(n):
            temp = commGen(ID,mod,pivot)
            if not temp == None:
                actives.append(ID)
                passives.append(temp)  
        mod = mod>>1
        pivot = pivot+1
        fnc(actives,passives,pivot,vargs)
        

            
def printlst(actives,passives,pivot,vargs): 
    '''
        The instantiation of the generalRun() virtual function for printing
        comm. data
    '''

    print "Round: "+str(pivot-2)
    realn = vargs
    for j in xrange(len(passives[0])):
        print ""
        print "\t",
        for i in xrange(len(actives)):
            if actives[i] < realn and passives[i][j]<realn:
                print str(actives[i]) + "->" +str(passives[i][j])+",",
    print ""
     
def printstd(n):
    '''
        This function prints the communication lists
        broken partitioned into rounds and comm. batches
    '''
    realn = n
    n = nextPow2(n)
    generalRun(n,printlst,realn) 



def drawConns(actives,passives,pivot,args):
    '''
        The instantiation of the generalRun() virtual function for drawing
        images
    '''
    [colors,draw,pts,wratio,res,allsteps,realn] = args
    draw.text(((wratio*res)*(pivot-2),1),"Round: "+str(pivot-2))
    for j in xrange(len(passives[0])):
        for i in xrange(len(actives)):
            t = (pts[actives[i]][0],pts[actives[i]][1])
            f= (pts[passives[i][j]][0],pts[passives[i][j]][1])
            if actives[i] < realn and passives[i][j]<realn:
                draw.line((t,f),fill=colors[(j) % len(colors)],width=2)
        if allsteps:
            for pt in xrange(len(pts)):
                pts[pt] = (pts[pt][0]+wratio*res,pts[pt][1])
                draw.text(pts[pt],str(pt))
                
    if j == 0: return
    
    if not allsteps:
        for pt in xrange(len(pts)):
            pts[pt] = (pts[pt][0]+wratio*res,pts[pt][1])
            draw.text(pts[pt],str(pt))
    
    


def drawImage(n,res=100,filename=None,allsteps=False):
    '''
        This function creates an image using the connection data
        it is written to conns.png by default, and has 300pixel resolution
        between points
    '''
    from math import log
    from PIL import Image
    import ImageDraw
    if filename==None:filename = str(n)+".png"
    realn = n
    n = nextPow2(n)
    wratio,hratio = mult(n)
    w = res*wratio
    h = res*hratio
    colors = ['red','sienna','green','cyan','blue','purple','yellow','turquoise','orange','brown','magenta','azure','khaki','chartreuse','coral','pink','tan','DarkOrchid','plum']

    
    if allsteps:img = Image.new('RGB', ((w)*(n-1),(h-res/2)))
    else:img = Image.new('RGB', ((w)*(int(log(n)/log(2))),(h-res/2)))
    draw = ImageDraw.Draw(img)
    
    pts = []
    for i in xrange(wratio):
        for j in xrange(hratio):
            pts.append((i*res,10+j*res))
    
    vargs = [colors,draw,pts,wratio,res,allsteps,realn]
    generalRun(n,drawConns,vargs)
    
    img.save(filename)


def roundChecker(actives,passives,pivot,vargs): 
    '''
        The instantiation of the generalRun() virtual function for checking
        comm. data
    '''
    alllsts = vargs
    
    for j in xrange(len(passives[0])):
        blocked = [0]*len(alllsts)
        
        for i in xrange(len(actives)):
            
            alllsts[actives[i]].append(passives[i][j])
            alllsts[passives[i][j]].append(actives[i])
            
            blocked[passives[i][j]] = blocked[passives[i][j]]+1
            blocked[actives[i]] = blocked[actives[i]]+1

        if sum(blocked)>len(blocked):
            for i in xrange(len(blocked)):
                if blocked[i]>1:print "In Round: "+str(pivot)+" cycle: "+str(j)+" node: " +str(i)+"was blocked"
     
def test(n):
    '''
        This function prints the communication lists
        broken partitioned into rounds and comm. batches
    '''
    
    n = nextPow2(n)
    alllsts = [[] for i in xrange(n)]
    generalRun(n,roundChecker,alllsts)  
    
    chk = sum(range(n))
    for i in range(n):
        if not sum(alllsts[i])==chk-i:print alllsts[i]

'''
    actual run stuff
'''
#printstd(32)


test(22)
drawImage(50,200)
#drawImage(64,300)
drawImage(22,50,allsteps=True)
#drawImage(16,100,allsteps=True)
#drawImage(8,100,allsteps=True)
#drawImage(4,100,allsteps=True)
#drawImage(8,100,allsteps=True)



