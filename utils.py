import numpy as np
from mystats import *

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
  
  
