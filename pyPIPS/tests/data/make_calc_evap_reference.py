"""Generate the reference output of DSDlib.calc_evap used by test_calc_evap_regression.py.

This was run ONCE against the original (pre-refactor) implementation of DSDlib.calc_evap
(commit cf256dd on branch evap_work, October 4, 2026) to capture its output. Do not re-run it
after the refactor unless you intend to change the reference behavior.

Usage: python make_calc_evap_reference.py
"""
from __future__ import annotations

import os

import numpy as np
from scipy.special import gammaln

import pyPIPS.DSDlib as dsd
import pyPIPS.thermolib as thermo


def make_inputs(n=400, seed=20261004):
    rng = np.random.default_rng(seed)
    T = rng.uniform(268., 308., n)  # K
    p = rng.uniform(65000., 103000., n)  # Pa
    RH = rng.uniform(20., 110., n)  # percent (includes supersaturated cases)
    mu = rng.uniform(-0.5, 15., n)
    lamda = 10.**rng.uniform(2.5, 4.0, n)  # m^-1
    Nt = 10.**rng.uniform(1., 4., n)  # m^-3
    N0 = Nt * np.exp((1. + mu) * np.log(lamda) - gammaln(1. + mu))  # m^-(4+mu)
    qv = thermo.calqv(np.minimum(RH, 100.) / 100., p, T)
    rho = thermo.calrho(p, thermo.caltheta(p, T), qv)
    return {'rho': rho, 'T': T, 'p': p, 'RH': RH, 'N0': N0, 'lamda': lamda, 'mu': mu}


if __name__ == '__main__':
    inputs = make_inputs()
    QREVP, COOL_RATE = dsd.calc_evap(inputs['rho'], inputs['T'], inputs['p'], inputs['RH'],
                                     inputs['N0'], inputs['lamda'], inputs['mu'])
    outfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'calc_evap_reference.npz')
    np.savez(outfile, QREVP=QREVP, COOL_RATE=COOL_RATE, **inputs)
    print(f"Wrote {outfile}")  # noqa: T201
