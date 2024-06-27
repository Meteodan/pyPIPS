from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotmodule as pm
from matplotlib.font_manager import FontProperties

font = {'size': 10}
mpl.rc('font', **font)

fontP = FontProperties()
fontP.set_size('small')


def empirical(zh, zdr):

    Nt = zh * 10.**(-0.0837 * (zdr**3.) + 0.702 * (zdr**2.) - 2.062 * zdr + 0.794)
    D0 = 0.0436 * (zdr**3.) - 0.216 * (zdr**2.) + 1.076 * zdr + 0.659
    W = zh * 10.**(-0.0493 * (zdr**3.) + 0.430 * (zdr**2.) - 1.542 * zdr - 3.019)
    R = zh * 10.**(-0.0363 * (zdr**3.) + 0.316 * (zdr**2.) - 1.178 * zdr - 1.964)

    return Nt, D0, W, R


def rad_emp_timeseries(dis, interp, obs, pstartindex, pstopindex, axparamdict,
                       DSDmidtimes, radmidtimes, image_dir, dis_name, name):

    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    fields = [dis[pstartindex:pstopindex + 1], interp[pstartindex:pstopindex + 1]]
    fieldparamdict1 = {
        'linestyle': '-',
        'color': 'g',
        'alpha': 0.5,
        'plotmin': 0,
        'label': r'rad obs'}
    fieldparamdict2 = {
        'linestyle': '-',
        'color': 'm',
        'alpha': 0.5,
        'plotmin': 0,
        'label': r'dis obs'}
    fieldparamdict3 = {
        'linestyle': '-',
        'color': 'c',
        'alpha': 0.5,
        'plotmin': 0,
        'label': r'rad interp'}
    fieldparamdicts = [fieldparamdict2, fieldparamdict3]
    xvals = [DSDmidtimes] * len(fields)
    ax1 = pm.plotmeteogram(ax1, xvals, fields, fieldparamdicts)
    axparamdicts = [axparamdict]
    ax1, = pm.set_meteogram_axes([ax1], axparamdicts)
    ax1.legend(loc='upper right', ncol=1, fancybox=True, shadow=False, prop=fontP)
    xvals.append(radmidtimes)
    fields.append(obs[pstartindex:pstopindex + 1])
    fieldparamdicts.append(fieldparamdict1)
    ax2 = pm.plotmeteogram(ax2, xvals, fields, fieldparamdicts)
    axes = [ax1, ax2]
    axparamdicts.append(axparamdict)
    ax2.legend(loc='upper left', ncol=2, fancybox=True, shadow=False, prop=fontP)
    axes = pm.set_meteogram_axes(axes, axparamdicts)
    plt.savefig(image_dir + 'meteograms/' + dis_name + '_' + name + '.png', dpi=300)
    plt.close(fig)


def retr_timeseries(obs, mm, retr_rad, retr_dis, pstartindex, pstopindex,
                    DSDmidtimes, axparamdict1, image_dir, dis_name, name):

    bias_dis = 100 * ((np.nansum(retr_dis - obs)) / np.nansum(obs))
    bias_rad = 100 * ((np.nansum(retr_rad - obs)) / np.nansum(obs))
    cc_dis = pd.DataFrame({'dis': retr_dis, 'obs': obs}).corr()
    cc_rad = pd.DataFrame({'rad': retr_rad, 'obs': obs}).corr()
    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(111)
    fields = [mm[pstartindex:pstopindex + 1], obs[pstartindex:pstopindex + 1],
              retr_rad[pstartindex:pstopindex + 1], retr_dis[pstartindex:pstopindex + 1]]
    fieldparamdict1 = {
        'linestyle': '-',
        'color': 'm',
        'alpha': 0.25,
        'plotmin': 0,
        'label': r'method of moments'}
    fieldparamdict2 = {
        'linestyle': '-',
        'color': 'k',
        'alpha': 0.25,
        'plotmin': 0,
        'label': r'observed'}
    fieldparamdict3 = {
        'linestyle': '-',
        'color': 'g',
        'alpha': 0.25,
        'plotmin': 0,
        'label': r'rad retrieved'}
    fieldparamdict4 = {
        'linestyle': '-',
        'color': 'c',
        'alpha': 0.25,
        'plotmin': 0,
        'label': r'dis retrieved'}
    fieldparamdicts = [fieldparamdict1, fieldparamdict2, fieldparamdict3, fieldparamdict4]
    xvals = [DSDmidtimes] * len(fields)
    ax1 = pm.plotmeteogram(ax1, xvals, fields, fieldparamdicts)
    axparamdicts = [axparamdict1]
    ax1, = pm.set_meteogram_axes([ax1], axparamdicts)
    if name in {'Nt', 'R'}:
        ax1.set_yscale('log')
    ax1.text(0.05, 0.93, f'Dis Retr. Bias ={bias_dis:2.2f}' + '%', transform=ax1.transAxes)
    ax1.text(0.05, 0.86, f'Rad Retr. Bias ={bias_rad:2.2f}' + '%', transform=ax1.transAxes)
    ax1.text(0.05, 0.79, f'Dis Retr. Corr Coeff ={cc_dis.iloc[0, 1]:2.3f}', transform=ax1.transAxes)
    ax1.text(0.05, 0.72, f'Rad Retr. Corr Coeff ={cc_rad.iloc[0, 1]:2.3f}', transform=ax1.transAxes)
    ax1.legend(
        bbox_to_anchor=(
            1.,
            1.),
        loc='upper right',
        ncol=1,
        fancybox=True,
        shadow=False,
        prop=fontP)
    plt.savefig(image_dir + 'meteograms/' + dis_name + '_' + name + '.png', dpi=300)
    plt.close(fig)


def dis_retr_timeseries(obs, retr_rad, retr_dis, pstartindex, pstopindex,
                        DSDmidtimes, axparamdict1, image_dir, dis_name, name):

    bias_dis = 100 * ((np.nansum(retr_dis - obs)) / np.nansum(obs))
    bias_rad = 100 * ((np.nansum(retr_rad - obs)) / np.nansum(obs))
    cc_dis = pd.DataFrame({'dis': retr_dis, 'obs': obs}).corr()
    cc_rad = pd.DataFrame({'rad': retr_rad, 'obs': obs}).corr()
    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(111)
    fields = [obs[pstartindex:pstopindex + 1],
              retr_rad[pstartindex:pstopindex + 1],
              retr_dis[pstartindex:pstopindex + 1]]
    fieldparamdict1 = {
        'linestyle': '-',
        'color': 'k',
        'alpha': 0.5,
        'plotmin': 0,
        'label': r'observed'}
    fieldparamdict2 = {
        'linestyle': '-',
        'color': 'g',
        'alpha': 0.5,
        'plotmin': 0,
        'label': r'rad retrieved'}
    fieldparamdict3 = {
        'linestyle': '-',
        'color': 'c',
        'alpha': 0.5,
        'plotmin': 0,
        'label': r'dis retrieved'}
    fieldparamdicts = [fieldparamdict1, fieldparamdict2, fieldparamdict3]
    xvals = [DSDmidtimes] * len(fields)
    ax1 = pm.plotmeteogram(ax1, xvals, fields, fieldparamdicts)
    axparamdicts = [axparamdict1]
    ax1, = pm.set_meteogram_axes([ax1], axparamdicts)
    ax1.text(0.05, 0.93, f'Dis Retr. Bias ={bias_dis:2.2f}' + '%', transform=ax1.transAxes)
    ax1.text(0.05, 0.86, f'Rad Retr. Bias ={bias_rad:2.2f}' + '%', transform=ax1.transAxes)
    ax1.text(0.05, 0.79, f'Dis Retr. Corr Coeff ={cc_dis.iloc[0, 1]:2.3f}', transform=ax1.transAxes)
    ax1.text(0.05, 0.72, f'Rad Retr. Corr Coeff ={cc_rad.iloc[0, 1]:2.3f}', transform=ax1.transAxes)
    ax1.legend(
        bbox_to_anchor=(
            1.,
            1.),
        loc='upper right',
        ncol=1,
        fancybox=True,
        shadow=False,
        prop=fontP)
    plt.savefig(image_dir + 'meteograms/' + dis_name + '_' + name + '.png', dpi=300)
    plt.close(fig)


def zh_zdr_timeseries(obs_rad, obs_dis, pstartindex, pstopindex,
                      DSDmidtimes, axparamdict1, image_dir, dis_name, name):

    bias = 100 * ((np.nansum(obs_rad - obs_dis)) / np.nansum(obs_dis))
    cc = pd.DataFrame({'dis': obs_dis, 'rad': obs_rad}).corr()
    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(111)
    fields = [obs_rad[pstartindex:pstopindex + 1], obs_dis[pstartindex:pstopindex + 1]]
    fieldparamdict1 = {
        'linestyle': '-',
        'color': 'g',
        'alpha': 0.5,
        'plotmin': 0,
        'label': r'radar obs'}
    fieldparamdict2 = {
        'linestyle': '-',
        'color': 'k',
        'alpha': 0.5,
        'plotmin': 0,
        'label': r'dis obs'}
    fieldparamdicts = [fieldparamdict1, fieldparamdict2]
    xvals = [DSDmidtimes] * len(fields)
    ax1 = pm.plotmeteogram(ax1, xvals, fields, fieldparamdicts)
    axparamdicts = [axparamdict1]
    ax1, = pm.set_meteogram_axes([ax1], axparamdicts)
    ax1.text(0.05, 0.93, f'Bias ={bias:2.2f}' + '%', transform=ax1.transAxes)
    ax1.text(0.05, 0.86, f'Corr Coeff ={cc.iloc[0, 1]:2.3f}', transform=ax1.transAxes)
    ax1.legend(
        bbox_to_anchor=(
            1.,
            1.),
        loc='upper right',
        ncol=1,
        fancybox=True,
        shadow=False,
        prop=fontP)
    plt.savefig(image_dir + 'meteograms/' + dis_name + '_' + name + '.png', dpi=300)
    plt.close(fig)


def one2one(obs, mm, retr_dis, retr_rad, image_dir, dis_name, name, obs_rain):
    if name == 'R':
        maxlim = 10**-1.
        minlim = 10**-4.
        yscale = 'log'
        label = 'R/Zh'
    if name == 'D0':
        maxlim = 4.0
        minlim = 0.0
        yscale = 'linear'
        label = 'D0'
    if name == 'Nt':
        maxlim = 100.0
        minlim = 10**-3.
        yscale = 'log'
        label = 'Nt/Zh'
    if name == 'W':
        maxlim = 10**-2.
        minlim = 10**-6.
        yscale = 'log'
        label = 'W/Zh'

    one_x = np.linspace(10**-8, 10**2)
    one_y = one_x
    bias_dis = 100 * ((np.nansum(retr_dis - obs)) / np.nansum(obs))
    bias_rad = 100 * ((np.nansum(retr_rad - obs)) / np.nansum(obs))
    cc_dis = pd.DataFrame({'dis': retr_dis, 'obs': obs}).corr()
    cc_rad = pd.DataFrame({'rad': retr_rad, 'obs': obs}).corr()
    fig1 = plt.figure(figsize=(8, 8))
    ax1 = fig1.add_subplot(111)
    plt.scatter(obs, retr_rad, c=obs_rain, marker='.', label='Rad Retrieval')
    plt.scatter(obs, retr_dis, c=obs_rain, marker='+', label='Dis Retreival')
    plt.scatter(obs, mm, c=obs_rain, label='Method Moments')
    cb = plt.colorbar()
    cb.set_label('Rainrate (mm/hr)')
    ax1.set_xlim(minlim, maxlim)
    ax1.set_xscale(yscale)
    ax1.set_ylim(minlim, maxlim)
    ax1.set_yscale(yscale)
    ax1.set_xlabel('Observed' + label)
    ax1.set_ylabel('Calculated' + label)
    ax1.plot(one_x, one_y, lw=2, color='k')
    ax1.text(0.6, 0.20, f'Dis Retr. Bias ={bias_dis:2.2f}' + '%', transform=ax1.transAxes)
    ax1.text(0.6, 0.15, f'Rad Retr. Bias ={bias_rad:2.2f}' + '%', transform=ax1.transAxes)
    ax1.text(0.6, 0.10, f'Dis Retr. Corr Coeff ={cc_dis.iloc[0, 1]:2.3f}', transform=ax1.transAxes)
    ax1.text(0.6, 0.05, f'Rad Retr. Corr Coeff ={cc_rad.iloc[0, 1]:2.3f}', transform=ax1.transAxes)
    plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=12.)
    plt.savefig(
        image_dir +
        'scattergrams/' +
        dis_name +
        '_one-to-' +
        name +
        '.png',
        dpi=200,
        bbox_inches='tight')
    plt.close(fig1)


def outer_one2one(obs, mm, retr_rad, image_dir, name, obs_rain):
    if name == 'R':
        maxlim = 10**-1.
        minlim = 10**-4.
        yscale = 'log'
        label = 'R/Zh'
    if name == 'D0':
        maxlim = 5.0
        minlim = 0.0
        yscale = 'linear'
        label = 'D0'
    if name == 'Nt':
        maxlim = 10**2.
        minlim = 10**-3.
        yscale = 'log'
        label = 'Nt/Zh'
    if name == 'W':
        maxlim = 10**-2.
        minlim = 10**-6.
        yscale = 'log'
        label = 'W/Zh'

    one_x = np.linspace(10**-8, 10**2)
    one_y = one_x
    bias_mm = 100 * ((np.nansum(mm - obs)) / np.nansum(obs))
    bias_retr = 100 * ((np.nansum(retr_rad - obs)) / np.nansum(obs))
    cc_mm = pd.DataFrame({'mm': mm, 'obs': obs}).corr()
    cc_retr = pd.DataFrame({'rad': retr_rad, 'obs': obs}).corr()
    fig1 = plt.figure(figsize=(8, 8))
    ax1 = fig1.add_subplot(111)
    plt.scatter(obs, mm, c=obs_rain, label='Method Moments')
    plt.scatter(obs, retr_rad, c=obs_rain, marker='.', label='Rad Retrieval')
    cb = plt.colorbar()
    cb.set_label('Rainrate (mm/hr)')
    ax1.set_xlim(minlim, maxlim)
    ax1.set_xscale(yscale)
    ax1.set_ylim(minlim, maxlim)
    ax1.set_yscale(yscale)
    ax1.set_xlabel('Observed' + label)
    ax1.set_ylabel('Calculated' + label)
    ax1.plot(one_x, one_y, lw=2, color='k')
    ax1.text(0.6, 0.20, f'Retr Bias ={bias_retr:2.2f}' + '%', transform=ax1.transAxes)
    ax1.text(0.6, 0.15, f'Moments Bias ={bias_mm:2.2f}' + '%', transform=ax1.transAxes)
    ax1.text(0.6, 0.10, f'Retr Corr Coeff ={cc_retr.iloc[0, 1]:2.3f}', transform=ax1.transAxes)
    ax1.text(0.6, 0.05, f'Moments Corr Coeff ={cc_mm.iloc[0, 1]:2.3f}', transform=ax1.transAxes)
    plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=12.)
    plt.savefig(image_dir + 'one-to-' + name + '.png', dpi=200, bbox_inches='tight')
    plt.close(fig1)


def scatters(obs, mm, retr_dis, retr_rad, ZDR_rad, ZDR, DSDmidtimes, image_dir, dis_name, name):  # noqa: ARG001

    if name == 'D0':
        ymin = 0.0
        ymax = 3.5
        ylabel = 'D0'
    if name == 'W':
        ymin = -6.0
        ymax = -1.0
        ylabel = 'log(W/Zh)'
    if name == 'R':
        ymin = -4.0
        ymax = -1.0
        ylabel = 'log(R/Zh)'
    if name == 'Nt':
        ymin = -2.5
        ymax = 1.5
        ylabel = 'log(Nt/Zh)'

    fig1 = plt.figure(figsize=(8, 8))
    ax1 = fig1.add_subplot(111)
    ax1.scatter(ZDR, retr_dis, color='c', marker='.', label='Dis Retrieval')
    ax1.scatter(ZDR, mm, color='m', marker='.', label='Method Moments')
    ax1.scatter(ZDR, obs, color='k', marker='.', label='Obs Disdrometer')
    ax1.scatter(ZDR_rad, retr_rad, color='g', marker='.', label='Rad Retrieval')
    ax1.set_xlim(-0.5, 3.5)
    ax1.set_ylim(ymin, ymax)
    ax1.set_xlabel('ZDR')
    ax1.set_ylabel(ylabel)
    plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=12.)
    plt.savefig(
        image_dir +
        'scattergrams/' +
        dis_name +
        '_' +
        name +
        '.png',
        dpi=200,
        bbox_inches='tight')
    plt.close(fig1)


def PIPS(obs_1A, obs_1B, obs_2A, obs_2B, ZDR_1A, ZDR_1B, ZDR_2A,
         ZDR_2B, ymin, ymax, image_dir, dis_name, name, ylabel):  # noqa: ARG001

    fig1 = plt.figure(figsize=(8, 8))
    ax1 = fig1.add_subplot(111)
    ax1.scatter(ZDR_1A, obs_1A, marker='.', label='PIPS 1A')
    ax1.scatter(ZDR_1B, obs_1B, marker='.', label='PIPS 1B')
    ax1.scatter(ZDR_2A, obs_2A, marker='.', label='PIPS 2A')
    ax1.scatter(ZDR_2B, obs_2B, marker='.', label='PIPS 2B')
    ax1.set_xlim(0.0, 4.0)
    ax1.set_ylim(ymin, ymax)
    ax1.set_xlabel('ZDR')
    ax1.set_ylabel(ylabel)
    plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=12.)
    plt.savefig(image_dir + 'scattergrams/' + name + '.png', dpi=200, bbox_inches='tight')
    plt.close(fig1)
