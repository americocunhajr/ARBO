import matplotlib
import numpy as np
from numpy import loadtxt

burnin = 100;
n_s = 7
var = 25000000
n_weeks = 52
dim = n_s + 1

# dataFile = "../inputs/datafile.txt"
dataFile = "sfp_qoi_seq.dat"
q = loadtxt(dataFile,comments="%")
q = q[burnin:]
# q = q.flatten()

sigma = np.sqrt(var)
# print(sigma)

filename = 'qoi-stats';
file = open(filename,'w')

#print(q.shape)
#print(range(q.shape[0]))
for i in range(n_weeks*dim):
    # (a, b) = divmod(i,n_s)
    # j = a*n_s + b
    f = [];
    # for j in range(len(q[:,0])):
    for j in range(q.shape[0]):
      mu = q[j,i];
      f.extend(np.random.normal(mu,sigma,100));
    # n, bins = np.histogram(q[:,i], bins = 200, normed = True)  #1.0*np.sqrt(len(q[:,i])))
    n, bins = np.histogram(f, bins = 200, density = True)#, normed = True)  #1.0*np.sqrt(len(q[:,i])))
    nn = np.cumsum(n*(bins[1]-bins[0]))
    q0 = max((bins[np.argmax(nn>.5)] + bins[np.argmax(nn>.5) +1])/2,0)
    q1 = max((bins[np.argmax(nn>.025)] + bins[np.argmax(nn>.025) +1])/2,0)
    q2 = max((bins[np.argmax(nn>.175)] + bins[np.argmax(nn>.25) +1])/2,0)
    q3 = max((bins[np.argmax(nn>.825)] + bins[np.argmax(nn>.75) +1])/2,0)
    q4 = max((bins[np.argmax(nn>.975)] + bins[np.argmax(nn>.975) +1])/2,0)
    # file.write(str(np.mean(q[:,i]))+' '+str(np.percentile(q[:,i],2.5))+' '+str(np.percentile(q[:,i],25))+' '+str(np.percentile(q[:,i],75))+' '+str(np.percentile(q[:,i],97.5))+'\n')
    #file.write(str(np.mean(q[:,i]))+' '+str(q1)+' '+str(q2)+' '+str(q3)+' '+str(q4)+'\n')
    file.write(str(q0)+' '+str(q1)+' '+str(q2)+' '+str(q3)+' '+str(q4)+'\n')

file.close()
