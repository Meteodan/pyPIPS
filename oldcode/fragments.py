# fragments.py
# Contains snippets of old code from different modules that may have a use somewhere at some point

# From pyPIPS.py:

if(pc.plot_radar and False):
    Zh_Cao = N.arange(20, 61, 1)
    Zdr_Cao = 10**((-2.6857 * 10**-4 * Zh_Cao**2) +
                   0.04892 * Zh_Cao - 1.4287)

    if(pc.plot_scat):
        if (not os.path.exists(ib.image_dir + 'scattergrams/')):
            os.makedirs(ib.image_dir + 'scattergrams/')
# 		fig1=plt.figure(figsize=(8,8))
# 		ax1=fig1.add_subplot(111)
# 		plt.title('Z vs. Zdr')
# 		ax1.scatter(dBZ, ZDR, color='k')
# 		ax1.scatter(radvars['dBZ'], radvars['ZDR'], color='r')
# 		ax1.plot(Zh_Cao,Zdr_Cao,lw=2)
# 		ax1.set_ylim(-2.0,6.0)
# 		ax1.set_xlim(20.0,60.0)
# 		ax1.set_xlabel('Zh in dBZ')
# 		ax1.set_ylabel('ZDR in dB')
# 		plt.savefig(image_dir+'scattergrams/'+dis_name+'scattergrams.png',dpi=200,bbox_inches='tight')
# 		plt.close(fig1)
#
# RELATIONS FROM CAO ET AL. 2008

    Zh_rad = pow(10., dBZ_rad / 10)

    # Radar measured
    Nt_rad_emp, D0_rad_emp, W_rad_emp, R_rad_emp = em.empirical(
        Zh_rad, ZDR_rad)

    # Disdrometer measured
    Nt_dis_emp, D0_dis_emp, W_dis_emp, R_dis_emp = em.empirical(Zh, ZDR)


# Timeseries and Figure 9a-c from Cao et al. and Figure 8a-c from Cao et al.
    name = 'R'
    axparamdict1 = {'majorxlocator': pc.locator,
                    'majorxformatter': pc.formatter,
                    'minorxlocator': pc.minorlocator,
                    'axeslimits': [[plotstarttime,
                                    plotstoptime],
                                   [0.0,
                                    250.0]],
                    'axeslabels': [pc.timelabel,
                                   r'Rain Rate']}
    em.retr_timeseries(
        PSD_df['intensity'].values,
        rainrate,
        R_rad_retr,
        R_dis_retr,
        pstartindex,
        pstopindex,
        PSDmidtimes,
        axparamdict1,
        ib.image_dir,
        dis_name,
        name)

    em.one2one(PSD_df['intensity'].values / Zh, rainrate / Zh, R_dis_retr / Zh,
               R_rad_retr / Zh_rad, ib.image_dir, dis_name, name)
    em.scatters(
        N.log10(
            PSD_df['intensity'].values / Zh),
        N.log10(
            rainrate / Zh),
        N.log10(
            R_dis_retr / Zh),
        N.log10(
            R_rad_retr / Zh_rad),
        ZDR_rad,
        ZDR,
        PSDmidtimes,
        ib.image_dir,
        dis_name,
        name)

    name = 'D0'
    axparamdict1 = {'majorxlocator': pc.locator,
                    'majorxformatter': pc.formatter,
                    'minorxlocator': pc.minorlocator,
                    'axeslimits': [[plotstarttime,
                                    plotstoptime],
                                   [0.0,
                                    5.0]],
                    'axeslabels': [pc.timelabel,
                                   r'D0']}
    em.retr_timeseries(
        D_med_disd,
        D_med_gam,
        D0_rad_retr,
        D0_dis_retr,
        pstartindex,
        pstopindex,
        PSDmidtimes,
        axparamdict1,
        ib.image_dir,
        dis_name,
        name)

    em.one2one(D_med_disd, D_med_gam, D0_dis_retr,
               D0_rad_retr, ib.image_dir, dis_name, name)
    em.scatters(D_med_disd, D_med_gam, N.array(D0_dis_retr), D0_rad_retr,
                ZDR_rad, ZDR, PSDmidtimes, ib.image_dir, dis_name, name)

    name = 'Nt'
    axparamdict1 = {'majorxlocator': pc.locator,
                    'majorxformatter': pc.formatter,
                    'minorxlocator': pc.minorlocator,
                    'axeslimits': [[plotstarttime,
                                    plotstoptime],
                                   [0.0,
                                    20000.0]],
                    'axeslabels': [pc.timelabel,
                                   r'Nt']}
    em.retr_timeseries(
        M0,
        Ntr_gam,
        Nt_rad_retr,
        Nt_dis_retr,
        pstartindex,
        pstopindex,
        PSDmidtimes,
        axparamdict1,
        ib.image_dir,
        dis_name,
        name)

    em.one2one(M0 / Zh, Ntr_gam / Zh, Nt_dis_retr / Zh,
               Nt_rad_retr / Zh_rad, ib.image_dir, dis_name, name)
    em.scatters(
        N.log10(
            M0 / Zh),
        N.log10(
            Ntr_gam / Zh),
        N.log10(
            Nt_dis_retr / Zh),
        N.log10(
            Nt_rad_retr / Zh_rad),
        ZDR_rad,
        ZDR,
        PSDmidtimes,
        ib.image_dir,
        dis_name,
        name)

    name = 'W'
    axparamdict1 = {'majorxlocator': pc.locator,
                    'majorxformatter': pc.formatter,
                    'minorxlocator': pc.minorlocator,
                    'axeslimits': [[plotstarttime,
                                    plotstoptime],
                                   [0.0,
                                    10.0]],
                    'axeslabels': [pc.timelabel,
                                   r'LWC']}
    em.retr_timeseries(
        LWC_disd,
        LWC_gam,
        W_rad_retr,
        W_dis_retr,
        pstartindex,
        pstopindex,
        PSDmidtimes,
        axparamdict1,
        ib.image_dir,
        dis_name,
        name)

    em.one2one(LWC_disd / Zh, LWC_gam / Zh, W_dis_retr / Zh,
               W_rad_retr / Zh_rad, ib.image_dir, dis_name, name)
    em.scatters(
        N.log10(
            LWC_disd / Zh),
        N.log10(
            LWC_gam / Zh),
        N.log10(
            W_dis_retr / Zh),
        N.log10(
            W_rad_retr / Zh_rad),
        ZDR_rad,
        ZDR,
        PSDmidtimes,
        ib.image_dir,
        dis_name,
        name)

    D0dict[dis_name + '_obs'] = D_med_disd
    ZDRdict[dis_name + '_obs'] = ZDR
    ZDRdict[dis_name + '_rad'] = Zh
    ZDRdict[dis_name] = dBZ
    Wdict[dis_name + '_obs'] = LWC_disd
    Rdict[dis_name + '_obs'] = PSD_df['intensity'].values
    Ntdict[dis_name + '_obs'] = M0

# name = 'D0'
# ymin = 0.0
# ymax = 5.0
# ylabel = 'D0'
# em.PIPS(D0dict['PIPS_1A_obs'],D0dict['PIPS_1B_obs'],D0dict['PIPS_2A_obs'],D0dict['PIPS_2B_obs'],ZDRdict['PIPS_1A_obs'],ZDRdict['PIPS_1B_obs'],ZDRdict['PIPS_2A_obs'],ZDRdict['PIPS_2B_obs'],ymin,ymax,image_dir,dis_name,name,ylabel)
#
# name = 'W'
# ymin = -6.0
# ymax = -1.0
# ylabel = 'log(W/Zh)'
# em.PIPS(N.log10(Wdict['PIPS_1A_obs']/ZDRdict['PIPS_1A_rad']),N.log10(Wdict['PIPS_1B_obs']/ZDRdict['PIPS_1B_rad']),N.log10(Wdict['PIPS_2A_obs']/ZDRdict['PIPS_2A_rad']),N.log10(Wdict['PIPS_2B_obs']/ZDRdict['PIPS_2B_rad']),ZDRdict['PIPS_1A_obs'],ZDRdict['PIPS_1B_obs'],ZDRdict['PIPS_2A_obs'],ZDRdict['PIPS_2B_obs'],ymin,ymax,image_dir,dis_name,name,ylabel)
#
# name = 'R'
# ymin = -5.0
# ymax = 0.0
# ylabel = 'log(R/Zh)'
# em.PIPS(N.log10(Rdict['PIPS_1A_obs']/ZDRdict['PIPS_1A_rad']),N.log10(Rdict['PIPS_1B_obs']/ZDRdict['PIPS_1B_rad']),N.log10(Rdict['PIPS_2A_obs']/ZDRdict['PIPS_2A_rad']),N.log10(Rdict['PIPS_2B_obs']/ZDRdict['PIPS_2B_rad']),ZDRdict['PIPS_1A_obs'],ZDRdict['PIPS_1B_obs'],ZDRdict['PIPS_2A_obs'],ZDRdict['PIPS_2B_obs'],ymin,ymax,image_dir,dis_name,name,ylabel)
#
# name = 'Nt'
# ymin = -4.0
# ymax = 2.0
# ylabel = 'log(Nt/Zh)'
# em.PIPS(N.log10(Ntdict['PIPS_1A_obs']/ZDRdict['PIPS_1A_rad']),N.log10(Ntdict['PIPS_1B_obs']/ZDRdict['PIPS_1B_rad']),N.log10(Ntdict['PIPS_2A_obs']/ZDRdict['PIPS_2A_rad']),N.log10(Ntdict['PIPS_2B_obs']/ZDRdict['PIPS_2B_rad']),ZDRdict['PIPS_1A_obs'],ZDRdict['PIPS_1B_obs'],ZDRdict['PIPS_2A_obs'],ZDRdict['PIPS_2B_obs'],ymin,ymax,image_dir,dis_name,name,ylabel)

# Figure 2 from Brandes et al. 2004
if(pc.plot_radar and False):
    one_x = N.linspace(0.0, 60.0)
    upper_y = N.exp(1.01 * 10**-4 * one_x**3 - 7.09 * 10**-
                    3 * one_x**2 + 2.38 * 10**-1 * one_x - 3.44)
    lower_y = N.exp(2.12 * 10**-4 * one_x**2 +
                    6.48 * 10**-2 * one_x - 3.87)

    fig1 = plt.figure(figsize=(8, 8))
    ax1 = fig1.add_subplot(111)
    ax1.scatter(ZDRdict['PIPS_1B'], ZDRdict['PIPS_1B_obs'],
                marker='.', label='PIPS 1B')
    ax1.scatter(ZDRdict['PIPS_2A'], ZDRdict['PIPS_2A_obs'],
                marker='.', label='PIPS 2A')
    ax1.scatter(ZDRdict['PIPS_2B'], ZDRdict['PIPS_2B_obs'],
                marker='.', label='PIPS 2B')
    ax1.set_xlim(0.0, 60.0)
    ax1.set_ylim(-1.0, 4.0)
    ax1.set_xlabel('Reflectivity (dBZ)')
    ax1.set_ylabel('ZDR (dB)')
    ax1.plot(one_x, upper_y, color='k')
    ax1.plot(one_x, lower_y, color='k')
    plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=8)
    plt.savefig(ib.image_dir + 'scattergrams/brandes.png',
                dpi=200, bbox_inches='tight')
    plt.close(fig1)

    # name = 'D0_obs'
    # ymin = 0.0
    # ymax = 5.0
    # ylabel = 'D0'
    # fig1=plt.figure(figsize=(8,8))
    # ax1=fig1.add_subplot(111)
    # ax1.scatter(ZDRdict['PIPS_1B_obs'], D0dict['PIPS_1Bobs'],color='m', marker='.', label='PIPS 1B')
    # ax1.scatter(ZDRdict['PIPS_2Aobs'], D0dict['PIPS_2Aobs'], color='k', marker='.', label='PIPS 2A')
    # ax1.scatter(ZDRdict['PIPS_2Bobs'], D0dict['PIPS_2Bobs'], color='g', marker='.', label='PIPS 2B')
    # ax1.set_xlim(0.0,4.0)
    # ax1.set_ylim(ymin,ymax)
    # ax1.set_xlabel('ZDR')
    # ax1.set_ylabel(ylabel)
    # plt.legend(loc='upper left',numpoints=1,ncol=1,fontsize=8)
    # plt.savefig(image_dir+'scattergrams/'+name+'.png',dpi=200,bbox_inches='tight')
    # plt.close(fig1)

    # Plot the lambda-mu relation and fit with 2nd order polynomial

    poly = N.polyfit(lamda, mu, 2)
    polynomial = N.poly1d(poly)

    xx = N.linspace(0.0, 30.0)
    yy = polynomial(xx)
    # yy2 = polynomial2(xx)

    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_subplot(111)
    plt.title('Shape-Slope Relation')
    ax1.scatter(lamda, mu, color='k')
    # ax1.scatter(Lam_retr,Mu_retr,color='g')
    ax1.plot(xx, yy, color='r')
    # ax1.plot(xx,yy2,color='b')
    # ax1.set_xlim(0.0,30.0)
    # ax1.set_ylim(-5.0,20.0)
    ax1.set_xlabel('Slope parameter')
    ax1.set_ylabel('Shape parameter')

    plt.savefig(ib.image_dir + 'scattergrams/' +
                'shape_slope.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    print(poly)
    # print(poly2)
