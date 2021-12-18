============================
mishchenko_brf python module
============================

Compute Bi-directional Reflection Function (BRF) in Python using Fortran Michael Mishchenko's radiative transfer solver.

The code computes the scalar bidirectional reflectance of a semi-infinite homogeneous slab composed of arbitrarily shaped, randomly oriented particles based on a rigorous numerical solution of the radiative transfer equation. 

The code is based on the publication: M. I. Mishchenko, J. M. Dlugach, E. G. Yanovitskij, and N. T. Zakharova, Bidirectional reflectance of flat, optically thick particulate layers: An efficient radiative transfer solution and applications to snow and soil surfaces, J. Quant. Spectrosc. Radiat. Transfer, 63, 409-432 (1999), https://doi.org/10.1016/S0022-4073(99)00028-X.

The Fortran code was downloaded from https://www.giss.nasa.gov/staff/mmishchenko/brf/ where some documentation is available.

License
-------

The Fortran code is a Free software for not-for-profit scientific research. The python binding is under MIT License.

Features
--------

* Compute spherical / diffuse albedo
* Compute direct albedo at 100 angles
* Compute BRDF at 100 x100 angles

Install
-------
    
    $ pip install git+https://github.com/ghislainp/mishchenko_brf

Example
-------

The main inputs are the single scattering albedo and legendre series coefficients of the phase function. The output is the brf, spherical albedo and direct albedo.


```python
import matplotlib.pyplot as plt
import numpy as np

from mishchenko_brf import brf

legendre = [1, 0, 0]  # isotropic phase function

single_scattering_albedo = 0.8

b = brf(single_scattering_albedo, legendre)

#
print(b.spectralalbedo())

#
theta = np.arange(0, 90)
alb = b.albedo(np.deg2rad(theta))

plt.figure()
plt.plot(theta, alb)

#
plt.figure()
plt.imshow(b.brf(np.deg2rad(theta), np.deg2rad(theta), phi=np.pi))

```


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage


It is an honor to make one of Michael Mishchenko's major contributions to snow optics available to a wider audience.
