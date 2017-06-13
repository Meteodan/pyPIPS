import matplotlib
import numpy as N
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.dates import *
from matplotlib.ticker import *
import glob
import modules.plotmodule as pm
from matplotlib.font_manager import FontProperties

font = {'size':10}
matplotlib.rc('font',**font)

fontP = FontProperties()
fontP.set_size('small')


def empirical(zh,zdr):

    Nt = zh * 10.**(-0.0837*(zdr**3.) + 0.702*(zdr**2.) - 2.062*zdr + 0.794)
    D0 = 0.0436*(zdr**3.) - 0.216*(zdr**2.) + 1.076*zdr + 0.659
    W = zh * 10.**(-0.0493*(zdr**3.) + 0.430*(zdr**2.) - 1.542*zdr - 3.019)
    R = zh * 10.**(-0.0363*(zdr**3.) + 0.316*(zdr**2.) - 1.178*zdr - 1.964)
    
    return Nt,D0,W,R
    
def rad_emp_timeseries(obs,emp,pstartindex,pstopindex,axparamdict,DSDmidtimes,radmidtimes,image_dir,dis_name,name):

    fig = plt.figure(figsize=(8,4))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    fields =[obs[pstartindex:pstopindex+1]]
    fieldparamdict1 = {'linestyle':'-','color':'k','alpha':0.5,'plotmin':0,'label':r'observed'}
    fieldparamdict2 = {'linestyle':'-','color':'r','alpha':0.5,'label':r'radar empirical'}
    fieldparamdicts = [fieldparamdict1]
    xvals = [DSDmidtimes]*len(fields)
    ax1 = pm.plotmeteogram(ax1,xvals,fields,fieldparamdicts)
    axparamdicts = [axparamdict]
    ax1, = pm.set_meteogram_axes([ax1],axparamdicts)
    ax1.legend(loc='upper right',ncol=1, fancybox=True, shadow=False, prop = fontP)
    xvals.append(radmidtimes)
    fields.append(emp[pstartindex:pstopindex+1])
    fieldparamdicts.append(fieldparamdict2)
    ax2 = pm.plotmeteogram(ax2,xvals,fields,fieldparamdicts)
    axes = [ax1,ax2]
    axparamdicts.append(axparamdict)
    ax2.legend(loc='upper left',ncol=1, fancybox=True, shadow=False, prop = fontP)
    axes = pm.set_meteogram_axes(axes,axparamdicts)
    plt.savefig(image_dir+'meteograms/'+dis_name+'_'+name+'.png',dpi=300)
    plt.close(fig)
    
def dis_emp_timeseries(obs,emp,mm,retr,pstartindex,pstopindex,DSDmidtimes,axparamdict1,image_dir,dis_name,name):

    fig = plt.figure(figsize=(8,4))
    ax1 = fig.add_subplot(111)
    fields = [emp[pstartindex:pstopindex+1],obs[pstartindex:pstopindex+1],mm[pstartindex:pstopindex+1],retr[pstartindex:pstopindex+1]]
    fieldparamdict1 = {'linestyle':'-','color':'m','alpha':0.5,'plotmin':0,'label':r'dis empirical'}
    fieldparamdict2 = {'linestyle':'-','color':'k','alpha':0.5,'plotmin':0,'label':r'observed'} 
    fieldparamdict3 = {'linestyle':'-','color':'c','alpha':0.5,'plotmin':0,'label':r'moments method'}
    fieldparamdict4 = {'linestyle':'-','color':'g','alpha':0.5,'plotmin':0,'label':r'dis retrieved'}
    fieldparamdicts = [fieldparamdict1,fieldparamdict2,fieldparamdict3,fieldparamdict4]
    xvals = [DSDmidtimes]*len(fields)
    ax1 = pm.plotmeteogram(ax1,xvals,fields,fieldparamdicts)
    axparamdicts = [axparamdict1]
    ax1, = pm.set_meteogram_axes([ax1],axparamdicts)
    ax1.legend(bbox_to_anchor=(1.,1.), loc='upper right',ncol=1, fancybox=True, shadow=False, prop = fontP)
    plt.savefig(image_dir+'meteograms/'+dis_name+'_'+name+'.png',dpi=300)
    plt.close(fig)

def dis_retr_timeseries(obs,emp,pstartindex,pstopindex,DSDmidtimes,axparamdict1,image_dir,dis_name,name):

    fig = plt.figure(figsize=(8,4))
    ax1 = fig.add_subplot(111)
    fields = [emp[pstartindex:pstopindex+1],obs[pstartindex:pstopindex+1]]
    fieldparamdict1 = {'linestyle':'-','color':'g','alpha':0.5,'plotmin':0,'label':r'dis retrieved'}
    fieldparamdict2 = {'linestyle':'-','color':'k','alpha':0.5,'plotmin':0,'label':r'observed'} 
    fieldparamdicts = [fieldparamdict1,fieldparamdict2]
    xvals = [DSDmidtimes]*len(fields)
    ax1 = pm.plotmeteogram(ax1,xvals,fields,fieldparamdicts)
    axparamdicts = [axparamdict1]
    ax1, = pm.set_meteogram_axes([ax1],axparamdicts)
    ax1.legend(bbox_to_anchor=(1.,1.), loc='upper right',ncol=1, fancybox=True, shadow=False, prop = fontP)
    plt.savefig(image_dir+'meteograms/'+dis_name+'_'+name+'.png',dpi=300)
    plt.close(fig)
    
def one2one(obs,emp,retr,mm,maxlim,minlim,image_dir,dis_name,name,yscale,label):

    one_x = N.linspace(10**-8,10**2)
    one_y = one_x

    fig1=plt.figure(figsize=(8,8))
    ax1=fig1.add_subplot(111)
    ax1.scatter(obs, emp, color='m', marker='.', label='Dis Empirical')
    ax1.scatter(obs, retr, color='g', marker='.', label='Dis Retreival')
    ax1.scatter(obs, mm, color='c', marker='.', label='Method Momemts')
    ax1.set_xlim(minlim,maxlim)
    ax1.set_xscale(yscale)
    ax1.set_ylim(minlim,maxlim)
    ax1.set_yscale(yscale)
    ax1.set_xlabel('Observed' + label)
    ax1.set_ylabel('Calculated' + label)
    ax1.plot(one_x,one_y,lw=2,color='k')
    plt.legend(loc='upper left',numpoints=1,ncol=1,fontsize=8)
    plt.savefig(image_dir+'scattergrams/'+dis_name+'_one-to-'+name+'.png',dpi=200,bbox_inches='tight')
    plt.close(fig1)
    
def scatters(obs,emp,retr,mm,ZDR,ymin,ymax,image_dir,dis_name,name,ylabel):

    bias_emp = 100 * (N.nansum(emp-obs)/N.nansum(obs))
    bias_retr = 100 * (N.nansum(retr-obs)/N.nansum(obs))
    fig1=plt.figure(figsize=(8,8))
    ax1=fig1.add_subplot(111)
    ax1.scatter(ZDR, emp, color='m', marker='.', label='Dis Empirical')
    ax1.scatter(ZDR, obs, color='k', marker='.', label='Obs Disdrometer')
    ax1.scatter(ZDR, retr, color='g', marker='.', label='Dis Retrieval')
    ax1.scatter(ZDR, mm, color='c', marker='.', label='Method Moments')
    ax1.set_xlim(0.0,4.0)
    ax1.set_ylim(ymin,ymax)
    ax1.set_xlabel('ZDR')
    ax1.set_ylabel(ylabel)
    ax1.text(0.8,0.95,'Emp. Bias =%2.2f'%bias_emp+'%',transform=ax1.transAxes)
    ax1.text(0.8,0.90,'Retr. Bias =%2.2f'%bias_retr+'%',transform=ax1.transAxes)
    plt.legend(loc='upper left',numpoints=1,ncol=1,fontsize=8)
    plt.savefig(image_dir+'scattergrams/'+dis_name+'_'+name+'.png',dpi=200,bbox_inches='tight')
    plt.close(fig1)