#Some handy functions
def extract_float(FileName,marker):
#the last occurance of the marker is extracted    
    import re
    #a = "average_pot =  "
    #def extract_float(FileName, marker):
    #extract a floating number, which follows marker
    #FileName = open('pos_dop_0/OUT.out')
    i=0
    for line in open(FileName):
        #marker = "average_pot =  "
        a = marker
        b = "[-+]?\d+\.\d+" 
        #print("a",a)
        #print("b",b)
        string = a+b
        match = re.findall(string, line)
        #print(match)
        if match != []:
           var = match[:]
           i=i+1
           #print(var)
    if i>0:       
       string = str(var)
       number = re.findall("[+-]?\d+\.\d+",string)
       number = number[0]
       number  = float(number)
       #type(number)
       return(number)
    else:
       return('artem: no matches found!')


def extract_e(FileName,marker):
#the last occurance of the marker is extracted    
    import re
    #a = "average_pot =  "
    #def extract_float(FileName, marker):
    #extract a floating number, which follows marker
    #FileName = open('pos_dop_0/OUT.out')
    i=0
    for line in open(FileName):
        #marker = "average_pot =  "
        a = marker
        b = "[-+]?\d+\.\d+[e][-+]?\d+" 
        #print("a",a)
        #print("b",b)
        string = a+b
        match = re.findall(string, line)
        #print(match)
        if match != []:
           var = match[:]
           i=i+1
           #print(var)
    if i>0:       
       string = str(var)
       number = re.findall("[-+]?\d+\.\d+[e][-+]?\d+",string)
       print('I am here', number[0])
       number = number[0]
       number  = float(number)
       #type(number)
       return(number)
    else:
       return('artem: no matches found!')

def extract_int(FileName,marker):
#the last occurance of the marker is extracted    
    import re
    i=0
    for line in open(FileName):
        a = marker
        b = "\d+"
        #print("a",a)
        #print("b",b)
        string = a+b
        match = re.findall(string, line)
        if match != []:
           var = match[:]
           i=i+1
    if i>0:       
       string = str(var)
       number = re.findall("\d+",string)
       number = number[0]
       number  = int(number)
       return(number)
    else:
       return('artem: no matches found!')

def average(array1d):
    #average defined as an average over last 50% of elements as a function of iteration number
    import numpy as np
    N = len(array1d)
    #N = 1000
    array2d = np.zeros(N)
    a = 0
    for i in range(0,N):
        a = array1d[int(i/2):i]
        array2d[i] = sum(a)/(len(a)+1)
    return array2d    

def time_average(t,f):
    #average defined as an average over last 50% of elements as a function of iteration number
    import numpy as np
    N = len(t)
    dt =  np.zeros(N)
    out = np.zeros(N)
    dt[0] = t[0]
    for i in range(1,N):
        dt[i] = t[i]-t[i-1]
    S = np.multiply(dt,f)
    #print(S)
    out[0]=f[0]
    t = np.array(t)
    for i in range(1,N):  
        j = (np.abs(0.5*t[i]-t)).argmin()#half-time idx 
        out[i] = sum(S[j:i])/(t[i]-t[j])
    return out  

def fun_bin1d_file(FileNameX, FileNameY, num_bins):
    from scipy.stats import binned_statistic
    import numpy as np

    #input
    #num_bins = 10000
    #FileNameX = 'sim_time.txt'
    #FileNameY = 'act_dop.txt'

    X = np.loadtxt(FileNameX,  dtype = float)
    Y = np.loadtxt(FileNameY,  dtype = float)

    statistic = binned_statistic(X, Y, statistic = 'mean',bins = num_bins)

    y = statistic[0]

    C = np.diff(statistic[1])
    x = statistic[1][1:]-C/2.0

    return x, y

#    plt.plot(x, y)
#    plt.savefig('artem.png')


def fun_bin1d(X, Y, num_bins):
    import numpy as np
    from scipy.stats import binned_statistic

    statistic = binned_statistic(X, Y, statistic = 'mean',bins = num_bins)

    y = statistic[0]

    C = np.diff(statistic[1])
    x = statistic[1][1:]-C/2.0

    return x, y

#    plt.plot(x, y)
#    plt.savefig('artem.png')

def fun_bin_timeav(X, Y, num_bins):
    import numpy as np
    from scipy.stats import binned_statistic

    num_i = len(X)
#    X0 = np.insert(X,0,0.0)#insert time 0 in the beggining
    dX = np.diff(X)#has a size of len(X)-1 dX[0] = X[1] - X[0]
    XY = np.zeros(num_i-1)
    XY = dX*Y[0:num_i-1]
    T = X[num_i - 1]#total time
    dT = T/num_bins
    XY = XY/dT

    statistic = binned_statistic(X[1:], XY, statistic = 'sum',bins = num_bins)

    y = statistic[0]

    C = np.diff(statistic[1])
    x = statistic[1][1:]-C/2.0

    return x, y

#    plt.plot(x, y)
#    plt.savefig('artem.png')

def fun_mid_time(A):
    import numpy as np
    a = A[len(A)-1]/2
    diff_array = np.abs(A-a)
    min_value = np.min(diff_array)
    idx = np.argwhere(min_value == diff_array)[0][0]

    return idx

def time_av_my(t,y,Nbin):

    import numpy as np
    import random

    Niter = len(t)
    dt = np.diff(t)

    T0 = min(t)#this is ueberfordert
    TN = max(t)#
    dTbin = (TN - T0)/Nbin

    IL  =   np.zeros(Nbin+1, dtype = np.int)
    IR  =   np.zeros(Nbin+1, dtype = np.int)
    T   =   np.linspace(T0,TN,Nbin+1)#bin borders

    for i in range(0, Nbin+1):
     dT = t - T[i]
     IL[i] = max(np.where(dT <= 0.0)[0])
     IR[i] = min(np.where(dT >= 0.0)[0])

    Y = np.zeros(Nbin, dtype = float)
    X = np.zeros(Nbin, dtype = float)

    dty = dt*y[0:Niter-1]
    for a in range(0, Nbin):
        b   =   np.sum(dty[ IR[a]:IL[a+1] ])
        c   =   ( t[IR[a]] - T[a]       ) * y[IL[a]]
        d   =   ( T[a+1]   - t[IL[a+1]] ) * y[IL[a+1]]
        Y[a]=   (b+c+d)/dTbin
        X[a]=   (T[a+1] + T[a])/2

    return X,Y


def get_slices_idx(A, n_slices):
        import numpy as np
        #return idxs of the sites, which x-coord lies in each of n_slices slices
        #1 - [ [x y z], [x, y, z], ... ]
        #2 - n slices to be cut
        #return: IDX; IDX[i] is the 1-d array of the corresponding idxs of the i-th slice
        Nidx = np.size(A,0)

        min_ = np.min(A[:,0])
        max_ = np.max(A[:,0])

        bins  =  np.linspace(min_,max_,n_slices)

        dx = (max_ - min_)/(bins - 1)
        bin_min_max = [bins[:] - 0.5, bins[:] + 0.5]

        IDX = np.array(np.zeros(n_slices), dtype = object)

        for i in np.arange(n_slices):
         IDX[i] = np.where(np.logical_and(A[:,0] >= bin_min_max[0][i], A[:,0] < bin_min_max[1][i]))[0]

        return(IDX)


def get_slices_idx_ext(A, n_slices, idx_of_interest):
        import numpy as np
        #return only idxs that belong to the idx_of_interest
        #return idxs of the sites, which x-coord lies in each of n_slices slices
        #1 - [ [x y z], [x, y, z], ... ]
        #2 - n slices to be cut
        #return: IDX; IDX[i] is the 1-d array of the corresponding idxs of the i-th slice
        Nidx = np.size(A,0)

        min_ = np.min(A[:,0])
        max_ = np.max(A[:,0])
        bins  =  np.linspace(min_,max_,n_slices)
        dx = (max_ - min_)/(bins - 1)#wtf?
        bin_min_max = [bins[:] - 0.5, bins[:] + 0.5]

        IDX = np.array(np.zeros(n_slices), dtype = object)

        for i in np.arange(n_slices):
         IDX[i] = idx_of_interest[np.where(np.logical_and(A[idx_of_interest,0] >= bin_min_max[0][i], A[idx_of_interest,0] < bin_min_max[1][i]))[0]]

        return(IDX)

def get_onset(E, DOS, a):
	import numpy as np
#onset of the distribution. Define the E_onset: \int_E[0]^E[i_onset](DOS) = int_-inf^inf DOS * a
#1 array of args
#2 array of fun values
# onset value a<~1
#DOS(E) vector to find an upper limit
	N		=	len(DOS)
	total_sum       =       np.sum(DOS)
	partial_sums    =       np.zeros(N, dtype = np.float)	
	for i in range(0, N):
        	partial_sums[i]         =       np.sum(DOS[0:i])
	partial_sums /= total_sum#normalization
	idx = np.argmin(np.abs(partial_sums - a))
#	print('who is closer to a = %f'%a, np.abs(partial_sums-a))
#	print('idx', idx)
#	print('E[idx]', E[idx])
	return E[idx]

