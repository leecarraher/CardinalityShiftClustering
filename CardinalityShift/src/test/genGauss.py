from numpy import random
part= 3
clu = 5
d = 8
print "x y z label"
for i in range(clu):
    variance = random.random()*4
    means = [random.random(1)*10 for b in range(d)] 
    
    for j in range(part):
        for k in xrange(d):
            pts = random.random()
            print (pts*variance+means[k])[0],
            
            
        print chr(i+97)
