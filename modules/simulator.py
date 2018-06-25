# simulator.py: a collection of functions related to the Parsivel simulator

import numpy as np
from scipy.stats import gamma, uniform
from . import disdrometer_module as dis


def samplegammaDSD(Nt, lamda, alpha, bins=None):
    """Randomly samples a gamma DSD given Nt, lamda, and alpha."""
    scale = 1. / lamda
    shape = alpha + 1.
    s = gamma.rvs(shape, scale=scale, size=int(Nt))
    # If desired, bin up the resulting samples into the given diameter ranges
    if(bins is None):
        return s
    else:
        ND_sample, _ = np.histogram(s, bins)
        return ND_sample / (bins[1:] - bins[:-1])


def create_random_gamma_DSD(Nt, lamda, alpha, Vt, sampling_length, sampling_width, Dl, Dmid, Dr,
                            Dmin=0, Dmax=20.0, sampling_interval=10., remove_margins=False,
                            verbose=False):
    """Given Nt, lamda, alpha, create a spatial distribution in a volume"""
    # First, determine the sampling volume. Use the sampling area A multiplied by the
    # depth that the fasted falling particle would cover in the time given by sampling_interval
    Dmax_index = np.searchsorted(Dr, Dmax / 1000.)
    if verbose:
        print "Dmax_index = ", Dmax_index
    Vtmax = Vt[Dmax_index]
    sampling_height = Vtmax * sampling_interval
    sampling_area = sampling_length * sampling_width
    sampling_volume = sampling_area * sampling_height  # Maximum sampling volume
    sampling_volumes_D = Vt[:Dmax_index + 1] * sampling_interval * \
        sampling_area  # Sampling volumes as a function of diameter D

    if verbose:
        print "sampling height = ", sampling_height
        print "sampling volume = ", sampling_volume

    # Next, create a uniform distribution of n=Nt*sampling_volume drops within
    # the volume given by sampling_volume
    n = int(Nt * sampling_volume)
    if verbose:
        print "number concentration = ", Nt
        print "number of particles in sampling volume = ", n
    xpos = uniform.rvs(0., sampling_length, ((n, 1)))
    ypos = uniform.rvs(0., sampling_width, ((n, 1)))
    zpos = uniform.rvs(0., sampling_height, ((n, 1)))

    # Next, determine the sizes of the drops by drawing the diameters randomly from a gamma
    # distribution given by n, lamda, and alpha

    diameters = samplegammaDSD(n, lamda, alpha)
    # Restrict diameters to be less than Dmax
    diameter_mask = np.where(diameters <= Dmax / 1000.)
    diameters = diameters[diameter_mask]
    xpos = xpos[diameter_mask]
    ypos = ypos[diameter_mask]
    zpos = zpos[diameter_mask]

    if verbose:
        print "number of particles less than Dmax = ", xpos.size

    # Now, figure out which drops in the volume won't fall through the sensor
    # area in the given time, and remove them

    velocities = dis.assignfallspeed(diameters * 1000.)
    depths = velocities * sampling_interval
    keepers = np.where(zpos.squeeze() - depths <= 0.)

    xpos = xpos[keepers]
    ypos = ypos[keepers]
    zpos = zpos[keepers]

    diameters = diameters[keepers]
    velocities = velocities[keepers]

    # Now, figure out which drops are margin fallers and flag them if desired
    # Below also masks out the margins on the ends of the beam, which may not be correct. One can
    # envision a drop falling and partially striking the top of the Parsivel, shearing it apart with
    # only part of the drop falling through the beam. OTOH, this should show up as a small drop
    # falling too fast, just like an "edge" margin faller, though...

    margins_xl = xpos.squeeze() - diameters / 2. < 0.
    margins_xr = xpos.squeeze() + diameters / 2. > sampling_length
    margins_yl = ypos.squeeze() - diameters / 2. < 0.
    margins_yr = ypos.squeeze() + diameters / 2. > sampling_width
    margin_mask = margins_xl | margins_xr | margins_yl | margins_yr

    if verbose:
        print "number of particles that fall through sampling area = ", xpos.size
        print "number of these that are margin fallers = ", margin_mask.sum()

    if remove_margins:
        if verbose:
            print "Removing margin fallers!"
        xpos = xpos[~margin_mask]
        ypos = ypos[~margin_mask]
        zpos = zpos[~margin_mask]
        diameters = diameters[~margin_mask]
        velocities = velocities[~margin_mask]

    # Bin up particles into the Parsivel diameter bins (later will add velocity too; right now all
    # particles are assumed to fall strictly along theoretical/empirical fall speed curve)
    Dedges = np.append(Dl, Dr[-1])
    Dedges = Dedges[:Dmax_index + 2]
    pcount_binned, _ = np.histogram(diameters, Dedges)
    # Compute ND of sample in original Parsivel diameter and velocity bins
    # Gets particle count/unit volume/diameter interval
    ND = pcount_binned / (sampling_volumes_D * (Dr[:Dmax_index + 1] - Dl[:Dmax_index + 1]))

    positions = np.hstack((xpos, ypos, zpos))
    sample_dict = {'positions': positions, 'diameters': diameters, 'velocities': velocities,
                   'ND': ND, 'margin_mask': margin_mask}
    return sample_dict