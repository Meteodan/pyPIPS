# pyPIPS_meteograms.py
#
# This script creates lookup tables for radar retrievals as a function of Z/ZDR over a certain
# range
import os
import argparse
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.ticker as ticker
import matplotlib.dates as dates
import matplotlib.pyplot as plt
import pyPIPS.parsivel_params as pp
import pyPIPS.radarmodule as radar
import pyPIPS.plotmodule as pm
import pyPIPS.pips_io as pipsio
import pyPIPS.utils as utils
import pyPIPS.PIPS as pips
import pyPIPS.parsivel_qc as pqc
import pyPIPS.timemodule as tm
import pyPIPS.DSDlib as dsd
import pyPIPS.polarimetric as dp

def DFcsv(output_dir, varname, var, dbz, zdr):
    var = np.array(var)
    var.shape = (len(dbz), len(zdr))
    df = pd.DataFrame(data=var, index=dbz, columns=zdr)
    df.to_csv(path_or_buf=output_dir + '/' + varname + '.csv', index_label='dBZ')

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

# Parse the command line options
description = "Computes lookup tables for radar retrievals"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('--coefficients', nargs=3, metavar=('c1', 'c2', 'c3'), type=float,
                    dest='coefficients',
                    help='coefficients of mu-lambda polynomial (in increasing order of exponent')
parser.add_argument('--ZH-range', nargs=3, metavar=('zb', 'ze', 'zi'), type=float, dest='ZH_range',
                    default=[0., 70., 0.1], help='beginning, ending, and interval for ZH range')
parser.add_argument('--ZDR-range', nargs=3, metavar=('zdb', 'zde', 'zdi'), type=float,
                    dest='ZDR_range', default=[0., 3., 0.01],
                    help='beginning, ending, and interval for ZDR range')
parser.add_argument('--scatt-file-path', dest='scatt_file_path',
                    help='path to scattering amplitude lookup table file')
parser.add_argument('--wavelength', type=float, dest='wavelength', default=10.7,
                    help='wavelength of radar in cm')
parser.add_argument('--output-dir', dest='output_dir', help='output directory')
args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

D_scatt, far_b, fbr_b, far_f, fbr_f = dp.readtmatrix(args.scatt_file_path)
fa2, fb2, fab, fba, far = dp.calbackscatterrain(far_b, fbr_b, far_f, fbr_f)

full_len = len(avg_diameter)
trunc_len = len(fa2)
D_trunc = avg_diameter[:trunc_len]
dD_trunc = bin_width[:trunc_len]

ZH = np.arange(args.ZH_range[0], args.ZH_range[1]+args.ZH_range[2], args.ZH_range[2])
ZDR = np.arange(args.ZDR_range[0], args.ZDR_range[1]+args.ZDR_range[2], args.ZDR_range[2])

RR = []
D0 = []
mu = []
lamda = []
N0 = []
Nt = []
W = []
sigma = []
Dm = []

for ZH_val in ZH:
    print("ZH = {:.2f}".format(ZH_val))
    RR2 = []
    D02 = []
    mu2 = []
    lamda2 = []
    N02 = []
    Nt2 = []
    W2 = []
    sigma2 = []
    Dm2 = []
    for ZDR_val in ZDR:
        retr_tuple = dsd.retrieval_Cao(ZH_val, ZDR_val, D_trunc, dD_trunc, fa2, fb2,
                                       args.wavelength, args.coefficients)

        RR2.append(retr_tuple[0])
        D02.append(retr_tuple[1])
        mu2.append(retr_tuple[2])
        lamda2.append(retr_tuple[3])
        N02.append(retr_tuple[4])
        Nt2.append(retr_tuple[5])
        W2.append(retr_tuple[6])
        sigma2.append(retr_tuple[7])
        Dm2.append(retr_tuple[8])

    RR.append(RR2)
    D0.append(D02)
    mu.append(mu2)
    lamda.append(lamda2)
    N0.append(N02)
    Nt.append(Nt2)
    W.append(W2)
    sigma.append(sigma2)
    Dm.append(Dm2)

DFcsv(args.output_dir, 'RR', RR, ZH, ZDR)
DFcsv(args.output_dir, 'D0', D0, ZH, ZDR)
DFcsv(args.output_dir, 'mu', mu, ZH, ZDR)
DFcsv(args.output_dir, 'lamda', lamda, ZH, ZDR)
DFcsv(args.output_dir, 'N0', N0, ZH, ZDR)
DFcsv(args.output_dir, 'Nt', Nt, ZH, ZDR)
DFcsv(args.output_dir, 'W', W, ZH, ZDR)
DFcsv(args.output_dir, 'sigma', sigma, ZH, ZDR)
DFcsv(args.output_dir, 'Dm', Dm, ZH, ZDR)
