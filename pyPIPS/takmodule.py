# takmodule.py -- variables and parameters used by the Takahashi bin scheme in COMMAS
import numpy as np

# Set up arrays for bin radius and volume for TAK scheme
rr_bin = 2.0e-6 * np.exp(np.arange(0., 34., 1.) / 4.329)
rbinwidth = 2. * rr_bin / 4.329
rs_bin = np.zeros((5, 21))
for i in range(5):
    rs_bin[i, :] = 2.0e-4 * np.exp(np.arange(0., 21., 1.) / 2.885)
# print rs_bin
hs_bin = 1.0e-4 * np.exp(np.tile(np.arange(0., 5., 1.), (21, 1)).transpose() /
                         (4. / np.log(2. * rs_bin / 2.0e-3)))    # Snow/Ice thickness
# print hs_bin
sbinwidth = 2. * rs_bin / 2.885
rg_bin = 2.0e-6 * np.exp(np.arange(0., 45., 1.) / 4.329)
rh_bin = rg_bin
gbinwidth = 2. * rg_bin / 4.329
hbinwidth = 2. * rh_bin / 4.329
vr_bin = (4.0 * np.pi / 3.0) * rr_bin**3.
ms_bin = (np.pi * (2.0 * rs_bin)**2. / 4.0) * \
    ((hs_bin - 1.0e-5) + 1.0e-5 * 100.)  # Mass of ice crystal
vg_bin = (4.0 * np.pi / 3.0) * rg_bin**3.
vh_bin = vg_bin
