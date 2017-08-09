# Creates CSV files for the radar retrievals. Would need to be updated anytime DSDretrieval_radar is changed
# command line: python pyPIPS_lookup.py

import matplotlib
import numpy as N
import pandas as pd
import modules.DSDretrieval_radar as DR
import modules.disdrometer_module as dis

### pulling in values from tmatrix stuff

filename = '/Users/bozell/pyPIPS/tmatrix/S-band/SCTT_RAIN_fw100.dat'
d,far_b,fbr_b,far_f,fbr_f = dis.readtmatrix(filename)
fa2,fb2,fab,fba,far = dis.calbackscatterrain(far_b,fbr_b,far_f,fbr_f)
min_diameter = dis.min_diameter
max_diameter = dis.max_diameter
intv = max_diameter-min_diameter
intv = intv[:N.size(fa2)]


dbz = N.arange(0.0,70.0,0.1)
zdr = N.arange(0.0,6.0,0.01)

def DFcsv(varname,var,dbz,zdr):
    var = N.array(var)
    var.shape = (len(dbz),len(zdr))
    df = pd.DataFrame(data = var, index=dbz, columns=zdr)
    df.to_csv(path_or_buf='lookups/'+varname+'.csv',index_label='dBZ')

R = []
D0 = []
MU = []
LAM = []
W = []
SIGM = []

for i in xrange(N.size(dbz)):
    R2 = []
    D0_2 = []
    MU2 = []
    LAM2 = []
    W2 = []
    SIGM2 = []
    for j in xrange(N.size(zdr)):
        R_retr,D0_retr,mu_retr,lam_retr,N0_retr,Nt_retr,W_retr,sigm_retr,Dm_retr,N_retr = DR.retrieve_DSD(dbz[i],zdr[j],d,fa2,fb2,intv)
        R2.append(R_retr)
        D0_2.append(D0_retr)
        MU2.append(mu_retr)
        LAM2.append(lam_retr)
        W2.append(W_retr)
        SIGM2.append(sigm_retr)
    R.append(R2)
    D0.append(D0_2)
    MU.append(MU2)
    LAM.append(LAM2)
    W.append(W2)
    SIGM.append(SIGM2)
    
DFcsv('R',R,dbz,zdr)
DFcsv('D0',D0,dbz,zdr)
DFcsv('mu',MU,dbz,zdr)
DFcsv('lam',LAM,dbz,zdr)
DFcsv('w',W,dbz,zdr)
DFcsv('sigm',SIGM,dbz,zdr)

