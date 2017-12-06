import pandas as pd
import numpy as np
import math
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import KFold


def fit_kde(data, cv=None):
    """
    Fit 1D density, representing number of points (Re/Im of visibility) with
    gaussian KDE.

    :param data:
        Numpy array.
    :return:
        Instance of ``sklearn.neighbors.KernelDensity`` class with best density
        estimate choosen by 5-fold CV.
    """

    params = {'bandwidth': np.logspace(-1.5, 0.5, 100)}
    grid = GridSearchCV(KernelDensity(), params, cv=cv)
    grid.fit(data[:, np.newaxis])

    return grid.best_estimator_


def k_s(alpha_0, m, n):
    return (1.0+alpha_0)*m + n - 3.0


def k_m(alpha_0, m, n):
    return ((3.0+2.0*alpha_0)*m+2.0*n-2.0)/(5.0+2.0*alpha_0)


def k_b(m):
    return 3.0*m-2.0


def alpha_s1(alpha_0, m, n):
    k_m_ = k_m(alpha_0, m, n)
    return (4.0+m-5.0*k_m_)/(2.0*k_m_)


def alpha_s2(alpha_0, m, n):
    k_s_ = k_s(alpha_0, m, n)
    k_b_ = k_b(m)
    return alpha_0 + k_s_/k_b_


def C1(alpha_0):
    # Charge of electron [C]
    e = 4.8 * 10 ** (-10)
    # Mass of electron [g]
    m = 9.109382 * 10 ** (-28)
    # Speed of light [cm / s]
    c = 3. * 10 ** 10
    return (np.sqrt(3.)*e**3./(m*c**2.*(2.*alpha_0+2.)))*(3.*e/(2.*np.pi*m*c))**alpha_0*math.gamma(0.5*alpha_0+11./6.)*math.gamma(0.5*alpha_0+1./6.)


def C2(alpha_0):
    # Charge of electron [C]
    e = 4.8 * 10 ** (-10)
    # Mass of electron [g]
    m = 9.109382 * 10 ** (-28)
    # Speed of light [cm / s]
    c = 3. * 10 ** 10
    return (m*c**2)**(2.*alpha_0)*(np.sqrt(3.)*e**3./(8.*np.pi*m))*(3.*e/(2.*np.pi*m**3*c**5))**(alpha_0+0.5)*math.gamma((6.*alpha_0+5.)/12.)*math.gamma((6.*alpha_0+25.)/12.)


def beta(Gamma):
    return math.sqrt(Gamma**2.-1.)/Gamma


def delta(Gamma, theta):
    return 1./(Gamma*(1.-beta(Gamma)*math.cos(theta)))


def comoving_transverse_distance(z, H_0=73.0, omega_M=0.3, omega_V=0.7,
                                 format="pc"):
    """
    Given redshift ``z``, Hubble constant ``H_0`` [km/s/Mpc] and
    density parameters ``omega_M`` and ``omega_V``, returns comoving transverse
    distance (see arXiv:astro-ph/9905116v4 formula 14). Angular diameter
    distance is factor (1 + z) lower and luminosity distance is the same factor
    higher.

    """
    from scipy.integrate import quad
    fmt_dict = {"cm": 9.26 * 10.0 ** 27.0, "pc": 3. * 10 ** 9, "Mpc": 3000.0,
                "Gpc": 3.0}

    result = (H_0 / 100.0) ** (-1.) * quad(lambda x: (omega_M * (1. + x ** 3) +
                                                      omega_V) ** (-0.5),
                                           0, z)[0]
    try:
        return fmt_dict[format] * result
    except KeyError:
        raise Exception('Format  \"pc\", \"cm\", \"Mpc\" or \"Gpc\"')


def s_obs(b1, n1, nu_obs, Gamma, fi, theta, z, m=1):
    """
    Observed flux (Jy) of BK jet for m={1, 2} and fixed n=2, alpha_0=0.5,
    alpha_si=alpha_s1 (that is nu_obs < nu_sM - characteristic frequency at the
    minimum distance of a core).

    :param b1:
        B at r=1pc [G].
    :param n1:
        N at r=1pc [cm^(-3)].
    :param nu_obs:
        Observed frequency [Hz].
    :param Gamma:
        Jet bulk motion Lorenz-factor.
    :param fi:
        Jet opening angle [rad].
    :param theta:
        LOS of jet.
    :param z:
        Redshift.
    :param m: (optional)
        Radial dependence of the magnetic field. (default: ``1``)
    """
    n = 2.0
    S1 = {1: -0.17, 2: -0.40}
    S2 = {1: 0.0, 2: -0.49}
    S3 = {1: 2.17, 2: 1.7}
    S4 = {1: 1.17, 2: 0.70}
    alpha_0 = 0.5
    delta_ = delta(Gamma, theta)
    beta_ = beta(Gamma)
    Dl9 = (1.+z)*comoving_transverse_distance(z, format="Gpc")
    number = 3.0*10**4*(6.2*10**18)**S1[m]*(6.9*10**7)**S2[m]
    return number*(C2(alpha_0))**S1[m]*Dl9**(-2.)*(1.+z)**(1.-alpha_0)*delta_**S3[m]*Gamma**(2.*S2[m])*beta_**(2.*S2[m])*fi**(2.+S1[m])*(1./math.sin(theta))**S1[m]*n1**(1.+S1[m])*b1**S4[m]*k_s(alpha_0, m, n)**(-1.)*nu_obs**(-alpha_s1(alpha_0, m, n))


# n_generate = 1000
# full_file = '/home/ilya/Dropbox/papers/accuracy/data/z_core_flux.txt'
# z_file = '/home/ilya/Dropbox/papers/accuracy/data/z.txt'
# names = ['source', 'epoch', 'class', 'z', 'S_core']
# df_z = pd.read_table(z_file, delim_whitespace=True, names=names,
#                      engine='python', usecols=[0, 1, 2, 3, 4])
# df = pd.read_table(full_file, delim_whitespace=True, names=names,
#                    engine='python', usecols=[0, 1, 2, 3, 4])
# zs = df_z['z'].values
# lzs = np.log(zs)
# # zs_mirrored = np.concatenate((zs, -zs), axis=0)
# cv = KFold(n_splits=5, shuffle=True)
# kde_lz = fit_kde(lzs, cv=cv)
#
# # Generate redshifts
# new_zs = np.exp(kde_lz.sample(n_generate)[:, 0])
#
# fluxes = df['S_core']
# lfluxes = np.log(fluxes)
# kde_lf = fit_kde(lfluxes, cv=cv)
#
# # Generate fluxes
# new_fluxes = np.exp(kde_lf.sample(n_generate)[:, 0])