#This is a library of functions to calculate various DSD-related parameters
# These were originally written in fortran but re-written using python and numpy

import numpy as N
from scipy.special import gamma as gamma_
from . import thermolib as thermo

def cal_Mp(p,lamda,Nt,alpha):
    """Given lamda, Nt, and alpha, computes the pth moment of the gamma distribution"""
    return (Nt/lamda**p)*(gamma_(1. + alpha + p)/gamma_(1. + alpha))

def cal_Dmpq(p,q,lamda,Nt,alpha):
    """Given two moments, p and q, and lamda,Nt, and alpha, compute the
       mean diameter associated with the ratio of p to q"""
    
    M_p = cal_Mp(p,lamda,Nt,alpha)
    M_q = cal_Mp(q,lamda,Nt,alpha)
    
    return (M_p/M_q)**(1./(p-q))

def cal_Nt(rhoa,q,N0,cx,alpha):
    """!   
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates number concentration at scalar points
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Dan Dawson 
    !  (02/06/2008) 
    !   
    !  MODIFICATION HISTORY:
    !  
    !  03/31/08 - converted intermediate calculations to double precision
    !             as well as a few of the input arguments.
    ! 
    !  04/24/09 - rewritten in Python/Numpy
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !     """

    #from numpy import *
    #from scipy.special import gamma as gamma_
    
    gamma1=gamma_(1.0+alpha)
    gamma4=gamma_(4.0+alpha)
    
    Ntx= (N0*gamma1)**(3.0/(4.0+alpha))* \
        ((gamma1/gamma4)*rhoa*q/cx)**((1.0+alpha)/(4.0+alpha))
    
    return Ntx

def cal_lamda(rhoa,q,Ntx,cx,alpha):
    """!
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates slope parameter lamda
    !
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Dan Dawson
    !  (02/06/2008)
    !
    !  MODIFICATION HISTORY:
    !  (03/31/2008)
    !  Converted intermediate calculations and arrays alpha and lamda to
    !  double precision. 
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !"""
    
    #from numpy import *
    #from scipy.special import gamma as gamma_
    
    gamma1 = gamma_(1.0+alpha)
    gamma4 = gamma_(4.0+alpha)

    lamda = ((gamma4/gamma1)*cx*Ntx/(rhoa*q))**(1.0/3.0)
    lamda = N.where(rhoa*q > 0.0,lamda,0.0)
    
    return lamda
    
def cal_N0(rhoa,q,Ntx,cx,alpha):
    """!   
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates intercept parameter and "effective" intercept parameter
    !            
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Dan Dawson 
    !  (02/06/2008) 
    !   
    !  MODIFICATION HISTORY:
    !   
    !  (03/26/2008)
    !  Recast N0 as a double precision variable, and used double precision for
    !  all intermediate calculations.  The calling subroutine should
    !  also define it as double precision.  For situations with large alpha,
    !  N0 can become very large, and loss of precision can result.
    !  Also tweaked the calculation of N0 a bit to avoid overflow, in keeping
    !  With Jason Milbrandt's calculation of N0 just before evaporation in 
    !  the multi-moment code.
    !
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !     """
    
    #from numpy import *
    #from scipy.special import gamma as gamma_
    
    gamma1 = gamma_(1.0+alpha)
    gamma4 = gamma_(4.0+alpha)
    
    lamda = cal_lamda(rhoa,q,Ntx,cx,alpha)
    lamda = lamda.astype(N.float64)
    try:
        alpha = alpha.astype(N.float64)
    except:
        pass
    
    N0 = Ntx*lamda**(0.50*(1.0+alpha))* \
        (1.0/gamma1)*lamda**(0.50*(1.0+alpha))
    
    N0_norm = N0*(((4.0+alpha)/lamda)** \
        alpha)*gamma4*(128.0/3.0)/ \
        ((4.0+alpha)**(4.0+alpha))
    N0 = N.where(lamda >= 0.0,N0,0.0)
    N0_norm = N.where(lamda >= 0.0,N0_norm,0.0)
    return N0,N0_norm

def cal_Dm(rhoa,q,Ntx,cx):
    """!   
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates mean-mass diameter
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Dan Dawson 
    !  (02/06/2008) 
    !   
    !  MODIFICATION HISTORY:
    !   
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !     """

    #from numpy import *

    Dm = (rhoa*q/(cx*Ntx))**(1./3.)
    Dm = N.where(Ntx > 0.0,Dm,0.0)
    
    return Dm

def cal_D0(rhoa,q,Ntx,cx,alpha):
    """!   
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates median volume diameter (m) for a gamma distribution
    !            Note, assumes spherical particles.
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Dan Dawson 
    !  (02/06/2008) 
    !   
    !  MODIFICATION HISTORY:
    !   
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !     """
    
    # Compute lamda
    lamda = cal_lamda(rhoa,q,Ntx,cx,alpha)
    
    return (3.672 + alpha)/lamda
    
def diag_alpha(rhoa,varid_qscalar,Dm):
    """!   
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates shape parameter alpha
    !-----------------------------------------------------------------------
    !     
    !  AUTHOR: Dan Dawson   
    !  (02/06/2008)         
    !   
    !  MODIFICATION HISTORY:
    !  (03/31/2008)
    !  Changed alpha array to double precision.
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !     """
    
    #from numpy import *
        
    alphaMAX = 80.0
    
    c1 = N.array([19.0,12.0,4.5,5.5,3.7])
    c2 = N.array([0.6,0.7,0.5,0.7,0.3])
    c3 = N.array([1.8,1.7,5.0,4.5,9.0])
    c4 = N.array([17.0,11.0,5.5,8.5,6.5])
    
    if(varid_qscalar == 'qr'):
        nq = 0
    elif(varid_qscalar == 'qi'): 
        nq = 1
    elif(varid_qscalar == 'qs'):
        nq = 2
    elif(varid_qscalar == 'qg'):
        nq = 3
    elif(varid_qscalar == 'qh'):
        nq = 4

    alpha = c1[nq]*N.tanh(c2[nq]*(1.e3*Dm-c3[nq]))+c4[nq]
    if(nq == 4):
        alpha=N.where(Dm > 0.008,1.e3*Dm-2.6,alpha)
        
    alpha = N.minimum(alpha,alphaMAX)
    
    return alpha

def solve_alpha(rhoa,cx,q,Ntx,Z):
    """!   
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates shape parameter alpha
    !-----------------------------------------------------------------------
    !     
    !  AUTHOR: Dan Dawson   
    !  (02/06/2008)         
    !   
    !  MODIFICATION HISTORY:
    !  (03/31/2008)
    !  Changed alpha array to double precision
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !     """
    
    #from numpy import *
    #from scipy.special import gamma as gamma_
    
    epsQ    = 1.e-14
    epsN    = 1.e-3
    epsZ    = 1.e-32
    
    alphaMax = 40.0
    
    tmp1= cx/(rhoa*q)
    g   = tmp1*Z*tmp1*Ntx
    g = N.where((q > epsQ) & (Ntx > epsN) & (Z > epsZ),g,-99.0)
    
    a = N.empty_like(q)
    a = N.where(g == -99.0,0.0,a)
    a = N.where(g >= 20.0,0.0,a)
    g2= g*g
    
    a = N.where((g < 20.0) & (g >= 13.31),3.3638e-3*g2 - 1.7152e-1*g + 2.0857e+0,a)
    a = N.where((g < 13.31) & (g >= 7.123),1.5900e-2*g2 - 4.8202e-1*g + 4.0108e+0,a)
    a = N.where((g < 7.123) & (g >= 4.200),1.0730e-1*g2 - 1.7481e+0*g + 8.4246e+0,a)
    a = N.where((g < 4.200) & (g >= 2.946),5.9070e-1*g2 - 5.7918e+0*g + 1.6919e+1,a)
    a = N.where((g < 2.946) & (g >= 1.793),4.3966e+0*g2 - 2.6659e+1*g + 4.5477e+1,a)
    a = N.where((g < 1.793) & (g >= 1.405),4.7552e+1*g2 - 1.7958e+2*g + 1.8126e+2,a)
    a = N.where((g < 1.405) & (g >= 1.230),3.0889e+2*g2 - 9.0854e+2*g + 6.8995e+2,a)
    a = N.where(g < 1.230,alphaMax,a)
    
    alpha = N.maximum(0.,N.minimum(a,alphaMax))
    
    return alpha

def calc_evap(rho,T,p,RH,N0,lamda,mu):
    """Computes bulk evaporation rate given the thermodynamic state of the air and gamma distribution
       parameters.  Much of this code was adapted from Jason Milbrandt's microphysics code.  Note,
       all thermodynamic variables are assumed to be in SI units."""

    # Thermodynamic constants
    
    Rd=287.0        #Gas constant for dry air (J/kg/K)
    Rv=461.51               #Gas constant for water vapor (J/kg/K)
    cp=1005.6       #Specific heat at constant pressure of dry air (J/kg/K)
    Lv=2.501e6      #Latent heat of vaporization (J/kg)
    Lf=3.34e5       #Latent heat of freezing (J/kg)
    Ls=Lv+Lf        #Latent heat of sublimation (J/kg)

    # Correction for fallspeeds at different pressures
    rho_ref = 1.225
    gamma_factor = N.sqrt(rho_ref/rho)
    #gamma_factor = 1.0

    cmr = (N.pi/6.)*1000.                                      # Constant in mass-diameter relation for rain

    # Calculate a bunch of quantities in the ventilation coefficient for the bulk evaporation rate
    Avx = 0.78                                              # Ventilation coefficients (Ferrier 1994)
    Bvx= 0.30 
    afr = 4854.0
    bfr = 1.0                                               # constants in fallspeed relation for rain from Ferrier 1994
    ffr = 195.0

    Cdiff = (2.2157e-5+0.0155e-5*(T-273.15))*1e5/p                  
    MUdyn = 1.72e-5*(393./(T+120.))*(T/273.16)**1.5     # Dynamic viscosity
    MUkin = MUdyn/rho                                       # Kinematic viscosity
    ScTHRD = (MUkin/Cdiff)**(1./3.)                         # Scorer parameter to the one-third power

    GR1 = gamma_(1.+mu)
    GR2 = gamma_(4.+mu)
    GR16 = gamma_(2.+mu)
    GR17 = gamma_(2.5+mu+0.5*bfr)

    cexr4 = 1.+mu
    cexr5 = 2.+mu
    cexr6 = 2.5+mu+0.5*bfr
    cexr9 = (GR2/GR1)*cmr

    # Now we can calculate the bulk ventilation coefficient for rain (Woohoo!)
    VENTr = Avx*GR16/(lamda**cexr5) + Bvx*ScTHRD*N.sqrt(gamma_factor*afr/MUkin)*GR17/(lamda+ffr)**cexr6

    #print 'VENTr = ',VENTr

    # Calculate the thermodynamic function in the denominator of the bulk evaporation rate

    Ka = 2.3971e-2 + 0.0078e-2*(T-273.15)                       # Thermal conductivity of air

    # Calculate saturation mixing ratio
    QSS = thermo.calqvs(p,T)

    #print 'QSS = ',QSS

    ABw = (Lv**2.)/(Ka*Rv*T**2.)+1./(rho*QSS*Cdiff)

    # Calculate saturation deficit S

    RH_frac = RH/100.0

    # Calculate water vapor mixing ratio
    qv = thermo.calqv(RH_frac,p,T)

    #print 'qv = ',qv

    S = qv/QSS-1.0

    #print 'S = ',S

    # With the above quantities, we can calculate the bulk evaporation rate (Yay!)

    QREVP = 2.*N.pi*S*N0*VENTr/ABw

    #print 'QREVP = ',QREVP

    # With QREVP we can calculate the latent cooling rate

    COOL_RATE = (Lv/cp)*QREVP
    
    return QREVP,COOL_RATE

def calc_VQR_ferrier(rhoa,q,Ntx,cx,alpha):
    """Calculates the mass-weighted rain fall speed using the formula in Ferrier (1994),
       assuming a gamma-diameter distribution and given q,N, and alpha"""
    
    afr= 4854.0  
    bfr= 1.0  
    ffr= 195.0 # Ferrier (1994)
    
    #First, compute slope parameter
    lamda = cal_lamda(rhoa,q,Ntx,cx,alpha)
    #Compute density correction factor
    gamfact = (1.204/rhoa)**0.5
    
    ckQr1 = afr*gamma_(4.0+alpha+bfr)/gamma_(4.0+ALFr)
    cexr1 = 4.0+alpha+bfr
    cexr2 = 4.00+alpha
    
    VQR = gamfact*ckQr1*lamda**(0.5*cexr2)/(lamda+ffr)**(0.5*cexr1)*lamda**(0.5*cexr2)/(lamda+ffr)**(0.5*cexr1)
    
    return VQR

def calc_VQR_gamvol(rhoa,q,Nt,cx,nu):
    """Calculates the mass-weighted rain fall speed for the NSSL scheme gamma-volume rain DSD"""
    
    #Compute density correction factor
    gamfact = (1.204/rhoa)**0.5
    #Compute mean volume
    vr = rhoa*q/(1000.*Nt)
    
    VQR = gamfact*(0.0911229*(1 + nu)**1.3333333333333333*gamma_(2. + nu) +              \
          5430.313059683277*(1 + nu)*vr**0.3333333333333333*                             \
          gamma_(2.333333333333333 + nu) -                                               \
          1.0732802065650471e6*(1 + nu)**0.6666666666666666*vr**0.6666666666666666*      \
          gamma_(2.6666666666666667 + nu) +                                              \
          8.584110982429507e7*(1 + nu)**0.3333333333333333*vr*gamma_(3 + nu) -           \
          2.3303765697228556e9*vr**1.3333333333333333*                                   \
          gamma_(3.333333333333333 + nu))/                                               \
          ((1 + nu)**2.333333333333333*gamma_(1 + nu))

    return VQR

def calc_VQR_gamdiam(rhoa,q,Nt,cx,alpha):
    """Calculates the mass-weighted rain fall speed for the NSSL scheme gamma-diameter rain DSD"""
    
    arx = 10.
    frx = 516.575 # raind fit parameters for arx*(1 - Exp(-fx*d)), where d is rain diameter in meters.
    
    #First, compute slope parameter
    lamda = cal_lamda(rhoa,q,Nt,cx,alpha)
    Dn = 1./lamda
    #Compute density correction factor
    gamfact = (1.204/rhoa)**0.5
    
    VQR = gamfact*arx*(1.0 - (1.0 + frx*Dn)**(-alpha - 4.0) )
    
    return VQR

def calc_VQG(rhoa,ax,bx,q,Ntx,cx,alpha):
    """Calculates mass-weighted graupel fall speed"""
    
    lamda = cal_lamda(rhoa,q,Ntx,cx,alpha)
    Dn = 1./lamda
    gamfact = (1.204/rhoa)**0.5
    VQG = gamfact*ax*(gamma_(4.0+alpha+bx)/gamma_(4.0+alpha))*Dn**0.5
    
    return VQG

def calc_VNG(rhoa,ax,bx,q,Ntx,cx,alpha):
    """Calculates number-weighted graupel fall speed"""
    
    lamda = cal_lamda(rhoa,q,Ntx,cx,alpha)
    Dn = 1./lamda
    gamfact = (1.204/rhoa)**0.5
    VNG = gamfact*ax*(gamma_(1.0+alpha+bx)/gamma_(1.0+alpha))*Dn**0.5
    
    return VNG


def power_mom(power,cx,t,q,moment):
    """!
    !-----------------------------------------------------------------------
    !
    ! PURPOSE:
    !
    ! Calculates moments of the PSD based on the Field et al. 2005 power law
    ! relations. Used for Thompson scheme.
    !
    ! AUTHOR:  Bryan Putnam, 4/16/2013
    !
    !-----------------------------------------------------------------------
    !"""

    T_c = t-273.16

    second_moment = (q/cx)

    if (power == 1):
       moment = second_moment
    else:
      log_a = 5.065339-.062659*T_c - 3.032362*power + 0.029469*T_c*power - \
              0.000285*(T_c**2) + 0.312550*(power**2) + 0.0000204*(T_c**2)*power + \
              0.003199*T_c*(power**2) + 0.000000*(T_c**3) - 0.015952*(power**3)

      a = 10.**log_a

      b = 0.476221 - 0.015896*T_c + 0.165977*power + 0.007468*T_c*power -   \
         0.000141*(T_c**2) + 0.060366*(power**2) + 0.000079*(T_c**2)*power + \
         0.000594*T_c*(power**2) + 0.000000*(T_c**3) - 0.003577*(power**3)

    moment = a*(second_moment)**b

    return moment
