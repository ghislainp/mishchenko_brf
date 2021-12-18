#!/usr/bin/env python

"""Tests for `mishchenko_brf` package."""

import numpy as np

from mishchenko_brf.lib.spher import mie


def test_mie():
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string

    ndistr = 3
    r1 = .245830
    r2 = 1.19417
    lam = 6300
    mrr = 1.530
    mri = 0.0080
    nk = 100
    n = 100
    np_ = 4
    ddelt = .1e-6
    npna = 0

    cext, csca, area, vol, al1 = mie(ndistr, r1, r2, lam, mrr, mri, n, np_, nk, npna, ddelt)

    expected_single_scattering_albedo, expected_legendre1 = results()
    print(al1.shape, expected_legendre1.shape)

    np.testing.assert_allclose(al1[:len(expected_legendre1)], expected_legendre1, atol=1e-5, rtol=0.)

    np.testing.assert_allclose(csca / cext, expected_single_scattering_albedo, atol=1e-6, rtol=0.)


def results():

    expected_single_scattering_albedo = 0.924351

    expected_AL1 = np.array(
        [
            1.00000,
            2.11107,
            2.80371,
            2.75088,
            2.78642,
            2.51483,
            2.37043,
            2.12491,
            1.94158,
            1.72589,
            1.54525,
            1.33856,
            1.18264,
            0.98046,
            0.85744,
            0.66452,
            0.57301,
            0.39838,
            0.33469,
            0.19091,
            0.15885,
            0.06616,
            0.06678,
            0.02026,
            0.02708,
            0.00370,
            0.00906,
            -0.00088,
            0.00188,
            0.00002,
            0.00001,
            0.00000,
            0.00000,
            0.00000,
        ]
    )
    return expected_single_scattering_albedo, expected_AL1
