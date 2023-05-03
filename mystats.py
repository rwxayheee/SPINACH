import numpy as np

# Variance of Shifted Data to avoid the catastrophic cancellation
def svar(mydata):
    # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    # Computing Shifted Data, let K (location parameter) be mean
    N = len(mydata)
    K = sum(mydata)/N
    
    sdata = [x - K for x in mydata]
    
    Sum1 = sum([x**2 for x in sdata])
    Sum2 = (sum(sdata))**2/N
    
    return (Sum1 - Sum2)/(N-1)

# Autocorrelation Function with default upper bound of lag = 200
def accor(mydata,lag_max=200):
    # https://scicoding.com/4-ways-of-calculating-autocorrelation-in-python/
    # modified definition of variance to 10.1021/ct0502864
    # computed as variance of shifted data, K = mean
    N = len(mydata)
    mymean = sum(mydata)/N
    myvar = svar(mydata)
    
    minlag = max(int(N/50),1)
    maxlag = min(N,int(lag_max))
    lags = range(minlag,maxlag)
    
    myacorr = {} 
    sdata = [x - mymean for x in mydata]
    
    for l in lags:
        c = 1 # Self correlation
        if (l > 0):
            tmp = [sdata[l:][i] * sdata[:-l][i] for i in range(len(mydata) - l)]
            c = sum(tmp)/(len(mydata) - l)/myvar
        myacorr[l] = c
        
    return myacorr

# Accumulated Average
def accavg(mydata):
    N = len(mydata)
    
    myacc=[]
    for i in range(1,N+1):
        myacc.append(sum(mydata[:i])/i)
    
    return myacc

# Utils

def en_parser(en_file):
    with open(en_file,"r") as f:
        line = f.readline()
        mydvdl = []
        while line:
            if 'L9' in line:
                mydvdl.append(line.split()[-1])
            line = f.readline()
    mydvdl = [float(x) for x in mydvdl[1:]]
    
    return mydvdl

def dvdl_stats(dvdldat):
    avgserial = accavg(dvdldat)
    return({"End-Point Average":avgserial[-1], "End-Point Variance":svar(avgserial)})

def make_integrand(Xs, en_list):
    fs = []
    vs = []
    for en in en_list:
        en_dvdl = en_parser(en)
        fs.append(dvdl_stats(en_dvdl)["End-Point Average"])
        vs.append(dvdl_stats(en_dvdl)["End-Point Variance"])
    return({"Xs":Xs, "fs":fs, "vs": vs})

def propvar(cubo, vs):
    return(sum(np.matmul(cubo.W_extra,np.array(vs).T)))

