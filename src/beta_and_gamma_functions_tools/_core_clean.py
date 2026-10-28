import math

# Lanczos approximation constants for g=7, n=9.
# g=7 is the smallest parameter set giving double-precision accuracy over the
# positive real axis, which is the only domain this library supports.
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


def log_gamma(x):
    """Natural log of the gamma function for real x > 0.

    The Lanczos approximation converges rapidly and reaches double-precision
    accuracy with a tiny coefficient table, which keeps the library dependency
    free.  Only positive real arguments are supported; the reflection formula
    is deliberately omitted to keep the scope small.
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
    return (
        0.5 * math.log(2.0 * math.pi)
        + math.log(s)
        + (z + 0.5) * math.log(t)
        - t
    )


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
