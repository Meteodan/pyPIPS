# DSDretrieval
# Adapted from code by Cao and Zhang

import numpy as np


def retrieve_DSD(Z, Zdr, d, fa2, fb2, intv, wavelength):
    # radar retrieval of DSD from Z(dBZ), and ZDR(dB) measurements
    # use constrained-gamma method

    na = len(Z)

    RAIN = []
    D_0 = []
    MU = []
    LAM = []
    N_0 = []
    NT = []
    LWC = []
    SIGM = []
    DM = []
    N_retr = []
    ND = []

    # initializing
    Dmx = 9     # maximum diameter of 9 mm

    epr = 80.205 + 1j * 17.167
    K2 = np.abs((epr - 1) / (epr + 2))**2.
#    K2 = 0.93       # dielectric constant
    wave = wavelength * 10.    # radar wavelength in mm
    v = -0.1021 + 4.932 * d - 0.9551 * d**2. + 0.07934 * d**3. - 0.002362 * \
        d**4.      # terminal velocity equation from Brandes et al. 2003
#    v = 3.778*d**0.67

# solve for lambda from the measured value of Zdr
    for n1 in range(0, na):     # running through each value in dBZ and ZDR array
        # iteration initialization
        dlm = 0.1
        mt = int(15 / dlm)        # these are the experimental lambda values to run through
        fp = -1.0       # used to test if observed Zdr greater than calculated
        lam1 = 1    # initial value
        zdr = 10.**(Zdr[n1] / 10)     # linear value
        refl = 10.**(Z[n1] / 10)      # linear value
        # dBZ = Z[n1]
        dbzdr = Zdr[n1]
        dmx = Dmx
        if np.all([refl >= 1, dmx >= 0.4, zdr >= 1]):     # when dBZ and ZDR are greater than 0
            if (dbzdr < 0.1):   # use empirical relations when ZDR < 0.1 dB from Cao et al. 2008
                lwc = refl * 10.**(-0.0411 * (dbzdr**3.) + 0.373 *
                                   (dbzdr**2.) - 1.439 * dbzdr - 2.9694)
                D0 = 0.0226 * (dbzdr**3.) - 0.0785 * (dbzdr**2.) + 0.8225 * dbzdr + 0.5811
                Dm = D0
                mum = np.nan
                lam = np.nan
                N0 = np.nan
                Nt = refl * 10.**(-0.0558 * (dbzdr**3.) + 0.5388 *
                                  (dbzdr**2.) - 1.8846 * dbzdr + 0.994)
                rain = refl * 10.**(-0.0307 * (dbzdr**3.) + 0.2768 *
                                    (dbzdr**2.) - 1.1138 * dbzdr - 1.9406)
                sigm = np.nan
#             elif (dbzdr > 3.0):     # ignore zdr values larger than 3dB
#                 rain = np.nan
#                 D0 = np.nan
#                 Dm = np.nan
#                 mum = np.nan
#                 lam = np.nan
#                 N0 = np.nan
#                 Nt = np.nan
#                 lwc = np.nan
#                 sigm = np.nan
            else:
                for n3 in range(1, mt):     # run through possible lambda values
                    tf = 0
                    lam = 0 + n3 * dlm      # assign number to lambda (0.1-30)
#                    mum = -0.0201*lam**2.+0.902*lam-1.718      # Cao et al. 2008
#                    mum = -0.016*lam**2. + 1.213*lam - 1.957    # Zhang et al 2001 relation
#                    mum =  -0.0264*lam**2. + 1.1073*lam-1.5206  # all IOP's
#                    mum  = -0.0324*lam**2. + 1.1931*lam - 1.7023 # all IOP's minus IOP_1A and
#                                                                 # plus IOP_1B_2 (with 2B)
#                    mum = -0.02014232*lam**2. + 1.05737156*lam - 2.09717828 ### FMCW Relation
#                    mum = -1.53068763 + 0.98197835*lam -0.01650507*lam**2. ### all FMCW, IOP data
                    # SATP relation, all FMCW and IOPs
                    mum = -1.0921 + 0.8172 * lam - 0.0126 * lam**2.
                    zh = 0.0
                    zv = 0.0
                    nd = d**mum * np.exp(-lam * d)
                    nd[d > Dmx] = 0.0       # if diameter larger than max diam, set  to 0.0
                    temp = fa2 * nd * intv
                    zh = np.sum(temp)
                    temp = fb2 * nd * intv
                    zv = np.sum(temp)
                    if(zv == 0.):
                        continue
                    rzdr = zh / zv        # calulated Zdr with this mu-lambda relation
                    fc = zdr - rzdr
                    # if the observed linear zdr is greater than the calculated one, exit loop with
                    # these values
                    if(fc * fp < 0):
                        lam = lam1 - np.abs(fp / (fc - fp)) * dlm  # not totally sure
                        tf = 1
                        break
                    else:
                        fp = fc
                        lam1 = lam
                # think this means a good lambda value was never found or that calculated zdr larger
                # than observed?
                if(tf == 0):
                    N0 = np.nan
                    mum = np.nan
                    lam = np.nan
                    D0 = np.nan
                    Dm = np.nan
                    sigm = np.nan
                    # these are relations taken from Brandes et al. 2004
                    Nt = 2.085 * refl * 10.**(0.728 * dbzdr**2. - 2.066 * dbzdr)
                    # intended for ZDR 0.3-3dB and calculated for
                    lwc = 5.589 * 10**-4. * refl * 10.**(0.223 * dbzdr**2. - 1.124 * dbzdr)
                    # lambda between 1-13 mm^-1
                    rain = 7.6 * 10**-3 * refl * 10.**(0.165 * dbzdr**2. - 0.897 * dbzdr)
                    if(dbzdr <= 0.35):
                        D0 = 0.717 + 1.479 * dbzdr - 0.725 * dbzdr**2. + 0.171 * dbzdr**3.
                        sigm = 0.163 + 0.519 * dbzdr - 0.0247 * dbzdr**2.
                else:
                    # SATP relation, all FMCW and IOPs
                    mum = -1.0921 + 0.8172 * lam - 0.0126 * lam**2.
                    zh0 = 0.0
                    # tkd = 0.0
                    Dm0 = 0.0
                    D3 = 0.0
                    D4 = 0.0
                    rain0 = 0.0

                    nd = d**mum * np.exp(-lam * d)
                    ND.append(nd)
                    nd[d > Dmx] = 0.0
                    zh0 = np.sum(fa2 * nd * intv)

                    D3 = np.sum(nd * d**3. * intv)
                    D4 = np.sum(nd * d**4. * intv)
                    rain0 = np.sum(np.pi * (0.6 * 10**-3) * v * d**3. * nd * intv)
                    Dm0 = np.sum(nd * intv)

                    Dm = D4 / D3
                    sigm2 = np.sum((d - Dm)**2. * nd * d**3. * intv)
                    sigm = (sigm2 / D3)**0.5

                    fhh2 = refl * np.pi**4. * K2 / 4.0 / wave**4.
                    N0 = fhh2 / zh0
                    Nt = Dm0 * N0
                    rain = rain0 * N0
                    lwc = (1 * 10**-3) * np.pi / 6 * N0 * D3

                    nd = N0 * d**mum * np.exp(-lam * d)
                    nd[d > Dmx] = 0.0
                    mtp1 = np.sum(nd * d**3. * intv)
                    mtp2 = 0.0
                    m = 0
                    mas = 0.0
                    if(mtp1 > 0):
                        while(mas < 0.5):
                            masp = mtp2 / mtp1
                            m = m + 1
                            mtp2 = mtp2 + nd[m] * d[m]**3. * intv[m]
                            mas = mtp2 / mtp1
                        D0 = d[m + 1] - (mas - 0.5) / (mas - masp) * intv[m + 1]
                    else:
                        D0 = 0.0
        else:
            rain = np.nan
            D0 = np.nan
            Dm = np.nan
            mum = np.nan
            lam = np.nan
            N0 = np.nan
            Nt = np.nan
            lwc = np.nan
            sigm = np.nan

        n_retr = N0 * d**mum * np.exp(-lam * d)
    #
        N_retr.append(n_retr)
        RAIN.append(rain)
        D_0.append(D0)
        MU.append(mum)
        LAM.append(lam)
        N_0.append(N0)
        NT.append(Nt)
        LWC.append(lwc)
        SIGM.append(sigm)
        DM.append(Dm)
    #

    return RAIN, D_0, MU, LAM, N_0, NT, LWC, SIGM, DM, N_retr
