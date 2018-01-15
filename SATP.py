import numpy as N
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl


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


radnpz_filename = 'SATP.npz'
radnpz_file = N.load(radnpz_filename)

Nc_bin_tarr = radnpz_file['Nc_bin']
Nc_bin_tarr = N.squeeze(Nc_bin_tarr,axis=0)
Nc_bin_tarr = N.swapaxes(Nc_bin_tarr,0,1)
R_tarr = N.squeeze(radnpz_file['R'],axis=0)
D0_tarr = N.squeeze(radnpz_file['D0'],axis=0)
mm = max(R_tarr)

rain_bins = []
D0_bins = N.arange(0.0,4.05,0.05)
D0_bins = N.array(D0_bins)
rainrate = 0.1
while rainrate < mm:
    rain_bins.append(rainrate)
    old_rainrate = rainrate
    rainrate = old_rainrate + (old_rainrate * 0.1)
R_bins = N.array(rain_bins)
print R_bins
print D0_bins
Nc_bin =[]
Nc_bin_avg = N.zeros((len(R_bins)-1,len(D0_bins)-1),dtype=object)
fig = plt.figure()
ax= fig.add_subplot(111)
for r in N.arange(len(R_bins)-1):
    for d in N.arange(len(D0_bins)-1):
        midpoint_R = N.median([R_bins[r],R_bins[r+1]])
        midpoint_D0 = N.median([D0_bins[d],D0_bins[d+1]])
#         print str(midpoint_R)+','+str(midpoint_D0) 
        Nc_bin = Nc_bin_tarr[N.logical_and(R_tarr>=R_bins[r],R_tarr<R_bins[r+1]) & N.logical_and(D0_tarr>=D0_bins[d],D0_tarr<D0_bins[d+1])]
        Nc_bin = N.where(Nc_bin == 0.0, N.nan, Nc_bin)
        Nc_bin_avg[r,d] = N.nanmean(Nc_bin,axis=0)
        if (len(Nc_bin) == 6.):
            ax.plot(avg_diameter,Nc_bin.T,color='0.7',alpha=0.5)
            ax.plot(avg_diameter,Nc_bin_avg[r,d],'k',label='(R,D0) = (%2.2f'%midpoint_R+',%2.2f'%midpoint_D0+')')
            ax.set_yscale('log')
            ax.set_ylim(10.**-1.0,10.**5.0)
            ax.set_ylabel('Number Concentration, # m^-3 mm^-1')
            ax.set_xlim(0.0,8.0)
            ax.set_xlabel('Diameter, mm')
            ax.tick_params(direction='in', length = 6, top = 'on',right = 'on')
            ax.legend(bbox_to_anchor=(1.,1.), loc='upper right',ncol=1, fancybox=True, shadow=False)
#plt.show()

fig = plt.figure()
ax = fig.add_subplot(211)
hist, xedges, yedges = N.histogram2d(R_tarr, D0_tarr, bins=(R_bins,D0_bins))
## add 2D color grid of bin sum
ax.set_xscale('log')
ax.pcolormesh(R_bins,D0_bins,hist.T)
# cax = fig.add_axes([0.90, 0.5, 0.02, 0.4])
# cb = mpl.colorbar.ColorbarBase(cax, spacing='proportional')

# Construct arrays for the anchor positions of the bars.
# Note: np.meshgrid gives arrays in (ny, nx) so we use 'F' to flatten xpos,
# ypos in column-major order. For numpy >= 1.7, we could instead call meshgrid
# with indexing='ij'.
ax = fig.add_subplot(212, projection='3d')
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
plt.show()



