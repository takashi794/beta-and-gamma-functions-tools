import math

# Lanczos approximation constants for g=7, n=9 (Numerical Recipes, commonly
# credited to Lanczos).  The choice g=7 is the smallest parameter set that
# delivers double-precision accuracy across the positive real axis, which is
# all this library supports.
_LANCZOS_G = 7.0
_LANCZOS_P = (
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
)

_PI = math.pi


def log_gamma(x):
    """Natural logarithm of the gamma function for real x > 0.

    The Lanczos approximation is used because it converges quickly and needs
    only a small table of coefficients to reach double-precision accuracy.

    Only positive real inputs are supported.  The reflection formula is not
    applied, so negative arguments and non-positive integers raise ValueError.
    The library is deliberately small rather than complete.
    """
    if not isinstance(x, (int, float)):
        raise TypeError("x must be a real number")
    if x <= 0:
        raise ValueError("x must be positive")
    z = x - 1.0
    s = _LANCZOS_P[0]
    for i in range(1, len(_LANCZOS_P)):
        s += _LANCZOS_P[i] / (z + i)
    t = z + _LANCZOS_G + 0.5
    return 0.5 * math.log(2.0 * _PI) + math.log(s) + (z + 0.5) * math.log(t) - t


def gamma(x):
    """Gamma function for real x > 0."""
    return math.exp(log_gamma(x))


def log_beta(a, b):
    """Natural log of the beta function B(a, b) for a > 0, b > 0."""
    if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
        raise TypeError("a and b must be real numbers")
    if a <= 0 or b <= 0:
        raise ValueError("a and b must be positive")
    return log_gamma(a) + log_gamma(b) - log_gamma(a + b)


def beta(a, b):
    """Beta function B(a, b) for a > 0, b > 0."""
    return math.exp(log_beta(a, b))
