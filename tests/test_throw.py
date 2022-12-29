import faulthandler
import unittest

import apts
import prrng

faulthandler.enable()


class Test_throw(unittest.TestCase):
    """
    Thrown particle
    """

    def test_simple(self):

        v0 = 100
        distribution = prrng.distribution.weibull
        parameters = [2]

        details = apts.ThrowParticleQuadratic(
            distribution_w=distribution,
            parameters_w=parameters,
            v0=v0,
        )

        wstop, vstop = apts.throw_particle_Quadratic(
            distribution_w=distribution,
            parameters_w=parameters,
            v0=[v0],
        )

        self.assertAlmostEqual(wstop, details.w[-1])
        self.assertAlmostEqual(vstop, details.v0[-1])


if __name__ == "__main__":

    unittest.main(verbosity=2)
