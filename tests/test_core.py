import math
import unittest

from beta_and_gamma_functions_tools import (
    log_gamma,
    gamma,
    log_beta,
    beta,
)


class TestLogGamma(unittest.TestCase):
    def test_integer_inputs(self):
        # gamma(n) = (n-1)!, so log_gamma(n) = log((n-1)!)
        self.assertAlmostEqual(log_gamma(1), math.log(math.factorial(0)), places=12)
        self.assertAlmostEqual(log_gamma(2), math.log(math.factorial(1)), places=12)
        self.assertAlmostEqual(log_gamma(5), math.log(math.factorial(4)), places=12)
        self.assertAlmostEqual(log_gamma(10), math.log(math.factorial(9)), places=12)

    def test_half_integer(self):
        # gamma(1/2) = sqrt(pi)
        self.assertAlmostEqual(log_gamma(0.5), 0.5 * math.log(math.pi), places=12)

    def test_known_value(self):
        # gamma(1.5) = sqrt(pi)/2;  log_gamma(1.5) = log(sqrt(pi)) - log(2)
        self.assertAlmostEqual(
            log_gamma(1.5), 0.5 * math.log(math.pi) - math.log(2.0), places=12
        )

    def test_large_argument(self):
        # Stirling-like growth: gamma(50) = 49!
        self.assertAlmostEqual(log_gamma(50), math.log(math.factorial(49)), places=8)

    def test_non_positive_raises(self):
        with self.assertRaises(ValueError):
            log_gamma(0)
        with self.assertRaises(ValueError):
            log_gamma(-1)


class TestGamma(unittest.TestCase):
    def test_gamma_half(self):
        self.assertAlmostEqual(gamma(0.5), math.sqrt(math.pi), places=12)

    def test_gamma_integer(self):
        self.assertAlmostEqual(gamma(5), math.factorial(4), places=8)

    def test_gamma_one(self):
        self.assertAlmostEqual(gamma(1), 1.0, places=14)


class TestLogBeta(unittest.TestCase):
    def test_symmetry(self):
        self.assertAlmostEqual(log_beta(2, 3), log_beta(3, 2), places=14)

    def test_known_value(self):
        # B(1, 1) = 1, log = 0
        self.assertAlmostEqual(log_beta(1, 1), 0.0, places=14)

    def test_half_values(self):
        # B(0.5, 0.5) = pi
        self.assertAlmostEqual(log_beta(0.5, 0.5), math.log(math.pi), places=12)

    def test_non_positive_raises(self):
        with self.assertRaises(ValueError):
            log_beta(0, 1)
        with self.assertRaises(ValueError):
            log_beta(1, -2)


class TestBeta(unittest.TestCase):
    def test_beta_identity(self):
        self.assertAlmostEqual(beta(1, 1), 1.0, places=14)

    def test_beta_known(self):
        # B(2, 3) = 1/12
        self.assertAlmostEqual(beta(2, 3), 1.0 / 12.0, places=12)

    def test_beta_half_half(self):
        self.assertAlmostEqual(beta(0.5, 0.5), math.pi, places=10)


if __name__ == "__main__":
    unittest.main()
