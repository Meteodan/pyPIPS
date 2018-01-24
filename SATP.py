import numpy as N
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import modules.disdrometer_module as dis
from numpy import ma as ma
import pyPIPScontrol as pc
from mpl_toolkits.axes_grid1 import make_axes_locatable

outer_image_dir='/Volumes/depot/dawson29/data/VORTEXSE/obsdata/jess_images/SATP/'

# Min diameter of bins (mm)
min_diameter_bins = [0.000,0.125,0.250,0.375,0.500,0.625,0.750,0.875,1.000,1.125,1.250,1.500,1.750,
                     2.000,2.250,2.500,3.000,3.500,4.000,4.500,5.000,6.000,7.000,8.000,9.000,10.000,
                     12.000,14.000,16.000,18.000,20.000,23.000]
max_diameter_bins = [0.125,0.250,0.375,0.500,0.625,0.750,0.875,1.000,1.125,1.250,1.500,1.750,2.000,
                     2.250,2.500,3.000,3.500,4.000,4.500,5.000,6.000,7.000,8.000,9.000,10.000,12.000,
                     14.000,16.000,18.000,20.000,23.000,26.000]
                     
# Average diameter of bins (mm)
avg_diameter_bins = [0.5*(x+y) for x,y in zip(min_diameter_bins,max_diameter_bins)]
avg_diameter = N.array(avg_diameter_bins)


radnpz_filename = 'SATP_all.npz'
radnpz_file = N.load(radnpz_filename)

Nc_bin_tarr = radnpz_file['Nc_bin']
# Nc_bin_tarr = N.squeeze(Nc_bin_tarr,axis=0)
# Nc_bin_tarr = N.swapaxes(Nc_bin_tarr,0,1)
# R_tarr = N.squeeze(radnpz_file['R'],axis=0)
# D0_tarr = N.squeeze(radnpz_file['D0'],axis=0)
R_tarr = radnpz_file['R']
D0_tarr = radnpz_file['D0']
mm = N.nanmax(R_tarr)
print mm
rain_bins = []
D0_bins = N.arange(0.0,4.05,0.05)
D0_bins = N.array(D0_bins)
rainrate = 0.1
while rainrate < mm:
    rain_bins.append(rainrate)
    old_rainrate = rainrate
    rainrate = old_rainrate + (old_rainrate * 0.1)
R_bins = N.array(rain_bins)

print len(R_bins) 
print len(D0_bins)
Nc_bin =[]
Nc_bin_avg = N.zeros((len(R_bins)-1,len(D0_bins)-1),dtype=object)

for r in N.arange(len(R_bins)-1):
    for d in N.arange(len(D0_bins)-1):
        midpoint_R = N.median([R_bins[r],R_bins[r+1]])
        midpoint_D0 = N.median([D0_bins[d],D0_bins[d+1]])
        Nc_bin = Nc_bin_tarr[N.logical_and(R_tarr>=R_bins[r],R_tarr<R_bins[r+1]) & N.logical_and(D0_tarr>=D0_bins[d],D0_tarr<D0_bins[d+1])]
        Nc_bin = N.where(Nc_bin == 0.0, N.nan, Nc_bin)
        Nc_bin_avg[r,d] = N.nanmean(Nc_bin,axis=0)
        if (len(Nc_bin) > 5):
            fig = plt.figure()
            ax= fig.add_subplot(111)
            ax.plot(avg_diameter,Nc_bin.T,color='0.7',alpha=0.5)
            ax.plot(avg_diameter,Nc_bin_avg[r,d],'k',label='(R,D0) = (%2.2f'%midpoint_R+',%2.2f'%midpoint_D0+')')
            ax.set_yscale('log')
            ax.set_ylim(10.**-1.0,10.**5.0)
            ax.set_ylabel('Number Concentration, # m^-3 mm^-1')
            ax.set_xlim(0.0,8.0)
            ax.set_xlabel('Diameter, mm')
            ax.tick_params(direction='in', length = 6, top = 'on',right = 'on')
            ax.legend(bbox_to_anchor=(1.,1.), loc='upper right',ncol=1, fancybox=True, shadow=False)
            plt.savefig(outer_image_dir+'DSD_R_'+str(midpoint_R)+'_D_'+str(midpoint_D0)+'.png',dpi=200,bbox_inches='tight')
            plt.close()


fig = plt.figure()
ax = fig.add_subplot(111)
hist, xedges, yedges = N.histogram2d(R_tarr, D0_tarr, bins=(R_bins,D0_bins))
## add 2D color grid of bin sum
ax.set_xscale('log')
ax.pcolormesh(R_bins,D0_bins,hist.T)
ax.set_xlabel('Rainrate')
ax.set_ylabel('D0')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cb = mpl.colorbar.ColorbarBase(cax, orientation = 'vertical')
plt.savefig(outer_image_dir+'2D_counts.png',dpi=200,bbox_inches='tight')
plt.close()

# Construct arrays for the anchor positions of the bars.
# Note: np.meshgrid gives arrays in (ny, nx) so we use 'F' to flatten xpos,
# ypos in column-major order. For numpy >= 1.7, we could instead call meshgrid
# with indexing='ij'.
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
xpos, ypos = N.meshgrid(xedges[:-1] + 0.25, yedges[:-1] + 0.25)
xpos = xpos.flatten('F')
ypos = ypos.flatten('F')
zpos = N.zeros_like(xpos)

# Construct arrays with the dimensions for the bars.
dx = 0.5 * N.ones_like(zpos)
dy = dx.copy()
dz = hist.flatten()
## add 3D bar graph, cannot get evenly spaced axis
ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')
ax.set_xlabel('Rainrate')
ax.set_ylabel('D0')
ax.set_zlabel('Nc_bin Count')
#plt.show()
plt.savefig(outer_image_dir+'3D_counts.png',dpi=200,bbox_inches='tight')
plt.close()

### fit with gamma distributions and calculate new mu-lambda relation
min_diameter = N.array(min_diameter_bins)
max_diameter = N.array(max_diameter_bins)
avg_size = avg_diameter
bin_width = max_diameter-min_diameter

mu = []
lam = []
for r in N.arange(len(R_bins)-1):
    for d in N.arange(len(D0_bins)-1):
        Nc_bin_use = Nc_bin_avg[r,d]
        Nc_bin_use = N.nan_to_num(Nc_bin_use)

        M2=N.sum(((avg_size[:]/1000.)**2.)*(1000.*Nc_bin_use)*bin_width[:]/1000.)
        M2 = ma.array(M2,dtype=N.float64)
        M4=N.sum(((avg_size[:]/1000.)**4.)*(1000.*Nc_bin_use)*bin_width[:]/1000.)
        M4 = ma.array(M4,dtype=N.float64)
        M6=N.sum(((avg_size[:]/1000.)**6.)*(1000.*Nc_bin_use)*bin_width[:]/1000.)
        M6 = ma.array(M6,dtype=N.float64)
        G =(M4**2.)/(M2*M6)
        G = ma.masked_invalid(G)
        mu_gam = ((7.-11.*G) - ((7.-11.*G)**2. - 4.*(G-1.)*(30.*G-12.))**(1./2.))/(2.*(G-1.))
        mu_gam = ma.masked_invalid(mu_gam)
        lamda_gam = ((M2*(mu_gam+3.)*(mu_gam+4.))/(M4))**(1./2.)
        lamda_gam = ma.masked_invalid(lamda_gam)
        mu.append(mu_gam)
        lam.append(lamda_gam/1000.)
        
## Plot the lambda-mu relation and fit with 2nd order polynomial 
lam = N.array(lam)
mu = N.array(mu)
lam = lam[~N.isnan(lam)]
mu = mu[~N.isnan(mu)]	
Lam = []
Mu = []
for n1 in xrange(0,len(lam)):
    lamda = lam[n1]
    if(lamda < 20.):
        Lam.append(lam[n1])
        Mu.append(mu[n1])
poly=N.polynomial.polynomial.polyfit(Lam,Mu,2)
polynomial=N.polynomial.polynomial.Polynomial(poly)
 
xx = N.linspace(0.0, 30.0)
yy = polynomial(xx)
y2 = -0.0201*xx**2. + 0.902*xx - 1.718
y3 = -0.016*xx**2. + 1.213*xx - 1.957
		
fig=plt.figure(figsize=(8,8))
ax1=fig.add_subplot(111)
plt.title('Shape-Slope Relation')
ax1.scatter(Lam,Mu, color='k', marker='.')
ax1.plot(xx,yy,label='Our Relation')
ax1.plot(xx,y2,label='Cao Relation')
ax1.plot(xx,y3,label='Zhang Relation')
ax1.set_xlabel('Slope parameter')
ax1.set_ylabel('Shape parameter')
ax1.text(0.05,0.90,'# of Points: %2.1f'%len(lam), transform=ax1.transAxes, fontsize=12.)
ax1.text(0.05,0.85,'%2.4f'%poly[2]+'*lam^2 + %2.4f'%poly[1]+'*lam + %2.4f'%poly[0], transform=ax1.transAxes, fontsize=12.)
plt.legend(loc='upper left',numpoints=1,ncol=3,fontsize=12.)
plt.savefig(outer_image_dir+'SATP_mu_lam.png',dpi=200,bbox_inches='tight')