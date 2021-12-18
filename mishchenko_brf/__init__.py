"""Top-level package for mishchenko_brf."""

__author__ = """Ghislain Picard"""
__email__ = 'ghislain.picard@univ-grenoble-alpes.fr'
__version__ = '0.9.0'

import os
import numpy as np
import scipy.interpolate
from warnings import warn


from mishchenko_brf.lib.refl import brf as fortran_brf
from mishchenko_brf.lib.spher import mie as fortran_mie


def brf(single_scattering_albedo, legendre_coefs, mmax=None, stdout=True):

    compiled_legendre_size = 700

    legendre_coefs = np.array(legendre_coefs)[:compiled_legendre_size].astype(np.float64)
    legendre_coefs = np.pad(legendre_coefs, (0, compiled_legendre_size - len(legendre_coefs)))

    if mmax is None:
        # same as the number of legendre_coefs
        warn("For computional time, it is recommended to set mmax to a reasonable value, not the default of 700.")
        mmax = len(legendre_coefs)

    fortran_brf_ = fortran_brf if stdout else mute_output(fortran_brf)

    nsep = 0  # could be one for a better numerical convergence but requires to do extra calculation (see interp.f)

    return BRDFResult(*fortran_brf_(single_scattering_albedo, mmax, legendre_coefs, nsep))


class BRDFResult(object):

    def __init__(self, ierr, mu, spherical_albedo, albedo, rzero, r):
        """create a BRDF results with the basic variables returned by refl.f.
"""
        self.ierr = ierr
        if ierr > 0:
            raise Exception(f"An error with code {ierr} occurred in refl fortran function")    
        self.mu = mu
        self._spherical_albedo = spherical_albedo
        self._albedo = albedo
        self.rzero = rzero
        self.r = r

    def spherical_albedo(self):
        """return the spherical/hemispheric albedo."""
        return self._spherical_albedo

    def albedo(self, theta):
        """return the albedo for the given incidence anglaes.

        :param theta: angle of incidence in radain
"""

        mu = np.cos(theta)
        return np.interp(mu, self.mu, self._albedo)

    def brf(self, theta_i, theta_v, phi, mmax=None):
        """compute the brf for given incidence and viewing angles
        :param theta_i: incidence angle in radian (float or 1d array)
        :param theta_v: viewing angle in radian (float or 1d array)
        :param phi: angle between incident and viewing angle (float or 1d array)
"""
        theta_i = np.atleast_1d(theta_i)
        theta_v = np.atleast_1d(theta_v)
        phi = np.atleast_1d(phi)

        assert len(theta_i.shape) <= 1
        assert len(theta_v.shape) <= 1
        assert len(phi.shape) <= 1

        def interp2d(z, mu_i, mu_v):
            interpolator = scipy.interpolate.interp2d(self.mu, self.mu, z, copy=False)
            return interpolator(mu_v, mu_i).reshape((len(mu_v), len(mu_i)))[:, :, np.newaxis]

        _mmax = self.r.shape[0]
        if mmax is not None:
            _mmax = min(mmax, _mmax)
        # Eq (2)

        def prepare_mu(theta):
            mu = np.cos(theta)
            order = mu.argsort()
            return mu[order], order

        mu_v, order_v = prepare_mu(theta_v)
        mu_i, order_i = prepare_mu(theta_i)

        phi = phi[np.newaxis, np.newaxis, :]

        res = interp2d(self.r[0], mu_i, mu_v) + \
             + 2 * sum(interp2d(self.r[m], mu_i, mu_v) * np.cos(m * phi) for m in range(1, _mmax))


        res = interp2d(self.r[0], mu_i, mu_v) + \
             + 2 * interp2d(self.r[1], mu_i, mu_v) * np.cos(1 * phi)

        # reorder the result
        res_out = np.empty((mu_v.size, mu_i.size, phi.size))
        res_out[order_v, :, :] = res
        res_out[:, order_i, :] = res_out.copy()

        return res_out.squeeze()

    # def brf_nearest_neighbor(self, theta_i, theta_v, phi, mmax=None):

    #     _mmax = self.r.shape[0]
    #     if mmax is not None:
    #         _mmax = min(mmax, _mmax)

    #     def prepare_mu(theta):
    #         mu = np.atleast_1d(np.cos(theta))
    #         order = mu.argsort()
    #         return mu[order], order

    #     mu_v, order_v = prepare_mu(theta_v)
    #     mu_i, order_i = prepare_mu(theta_i)

    #     ii = np.searchsorted(self.mu, mu_i)
    #     iv = np.searchsorted(self.mu, mu_v)

    #     ii[ii == len(self.mu)] = len(self.mu) - 1
    #     iv[iv == len(self.mu)] = len(self.mu) - 1

    #     res = self.r[0, iv, ii] + 2 * sum(self.r[m, iv, ii] * np.cos(m * phi) for m in range(1, _mmax))

    #     NotImplementedError("add the 'reorder' part. See above")
        #     if len(res.shape) == 2:
        #     res_out[order_v, order_i] = res.squeeze()
        # elif len(res.shape) == 1:
        #     if len(order_i) > 1:
        #         res_out[order_i] = res.squeeze()
        #     elif len(order_v) > 1:
        #         res_out[order_v] = res.squeeze()
        #     else:
        #         res_out = res
        # else:
        #     raise NotImplementedError("non implemented")

    #     return res


def mie(lam, r1, mrr, mri, ndistr=1, r2=None, n=10, np=4, nk=10, npna=0, ddelt=1e-5):
    """

    :param  ndistr: specifies the type of particle size distribution as follows:

     NDISTR = 1 - modified gamma distribution
          [Eq. (5.242) of Ref. 1]
              AA=alpha
              BB=r_c
              GAM=gamma
     NDISTR = 2 - log normal distribution
          [Eq. (5.243) of Ref. 1]
              AA=r_g
              BB=[ln(sigma_g)]**2
     NDISTR = 3 - power law distribution
          [Eq. (5.244) of Ref. 1]
               AA=r_eff (effective radius)
               BB=v_eff (effective variance)
               Parameters R1 and R2 (see below) are calculated
               automatically for given AA and BB
     NDISTR = 4 - gamma distribution
          [Eq. (5.245) of Ref. 1]
               AA=a
               BB=b
     NDISTR = 5 - modified power law distribution
          [Eq. (5.246) of Ref. 1]
               BB=alpha
     NDISTR = 6 - bimodal volume log normal distribution
              [Eq. (5.247) of Ref. 1]
              AA1=r_g1
              BB1=[ln(sigma_g1)]**2
              AA2=r_g2
              BB2=[ln(sigma_g2)]**2
              GAM=gamma

   :param R1, R2: minimum and maximum
         radii in the size distribution for NDISTR=1-4 and 6.
         R1 and R2 are calculated automatically
         for the power law distribution with given r_eff and v_eff
         but must be specified for other distributions.
         For the modified power law distribution (NDISTR=5), the
         minimum radius is 0, R2 is the maximum radius,
         and R1 is the intermediate radius at which the
         n(r)=const dependence is replaced by the power law
         dependence.

   :param LAM: wavelength of the incident light in the surrounding medium

   - Important:  r_c, r_g, r_eff, a, LAM, R1, and R2
                 must be in the same units of length (e.g., microns)

   :param MRR, MRI: real and imaginary parts of the relative refractive
                   index (MRI must be non-negative)

   :param N: number of integration subintervals on the interval (R1, R2)
         of particle radii
   :param NP: number of integration subintervals on the interval (0, R1)
         for the modified power law distribution
   :param NK: number of Gaussian division points on each of the
          integration subintervals

   :param NPNA: number of scattering angles at which the scattering
        matrix is computed
        (see the PARAMETER statement in subroutine MATR).
        The corresponding scattering angles are given by
        180*(I-1)/(NPNA-1) (degrees), where I numbers the angles.
        This way of selecting scattering angles can be easily changed
        at the beginning of subroutine MATR by properly modifying
        the following lines:

          N=NPNA
          DN=1D0/DFLOAT(N-1)
          DA=DACOS(-1D0)*DN
          DB=180D0*DN
          TB=-DB
          TAA=-DA
          DO 500 I1=1,N
             TAA=TAA+DA
             TB=TB+DB

        and leaving the rest of the subroutine intact.
        This flexibility is provided by the fact that
        after the expansion coefficients AL1,...,BET2 (see below)
        are computed by subroutine SPHER, the scattering matrix
        can be computed for any set of scattering angles without
        repeating Lorenz-Mie calculations.

   :param DDELT: desired numerical accuracy of computing the scattering
        matrix elements

  :returns:
      - IERR - error code (0 =no error)
      - CEXT and CSCA - average extinction and scattering cross sections
             per particle

      - ALBEDO - single-scattering albedo
      - AREA - average projected area per particle
      - VOL  - average volume per particle
      - ALPHA1 - coefficients appearing in the expansions
             of the elements of the normalized scattering matrix in
             generalized spherical functions [Eqs. (4.75)-(4.80) of Ref. 1].
"""
    if r2 is None:
        r2 = r1 * 1.0001
    return MieResult(*fortran_mie(int(ndistr), float(r1), float(r2), float(lam),
                                  float(mrr), float(mri),
                                  int(n), int(np), int(nk), int(npna), float(ddelt)))


class MieResult(object):

    def __init__(self, ierr, cext, csca, area, vol, al1):

        self.ierr = ierr
        if ierr > 0:
            raise Exception(f"An error with code {ierr} occurred in mie fortran function")
        self.cext = cext
        self.csca = csca
        self.single_scattering_albedo = csca / cext if cext > 0 else np.nan
        self.area = area
        self.vol = vol
        self.legendre1 = al1
        self.g = al1[1] / 3


def mute_output(func):

    def inner(*args, **kwargs):

        null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # save the current file descriptors to a tuple
        save = os.dup(1), os.dup(2)
        # put /dev/null fds on 1 and 2
        os.dup2(null_fds[0], 1)
        os.dup2(null_fds[1], 2)

        # *** run the function ***
        res = func(*args, **kwargs)
        # restore file descriptors so I can print the results
        os.dup2(save[0], 1)
        os.dup2(save[1], 2)
        # close the temporary fds
        os.close(null_fds[0])
        os.close(null_fds[1])

        return res

    return inner
