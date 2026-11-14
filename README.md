# beta-and-gamma-functions-tools

Computes the log-gamma and beta functions for positive real arguments using the Lanczos approximation. No third-party dependencies; standard-library `math` only.

```python
from beta_and_gamma_functions_tools import log_gamma, gamma, log_beta, beta

print(log_gamma(5.0))   # ~3.178 (i.e. ln 24)
print(gamma(0.5))      # ~1.7725 (i.e. sqrt pi)
print(log_beta(2, 3))  # ~-2.485 (i.e. ln(1/12))
print(beta(2, 3))      # ~0.0833 (i.e. 1/12)
```

## Why this exists

A need for a small, self-contained log-gamma / beta routine that runs anywhere the Python standard library runs — no `scipy`, no `mpmath`, no build step. The Lanczos approximation with `g=7` and nine coefficients is chosen because it delivers double-precision accuracy on the positive real axis with a tiny constant table and no iterative refinement.

## Edge cases

Only **positive real** arguments are supported. The reflection formula for negative arguments is deliberately not implemented, so `log_gamma(-0.5)` raises `ValueError` rather than returning a wrong or surprising value. Non-positive integers — the poles of gamma — also raise. Complex inputs are out of scope; if you need them, use a library that implements them.

## API

- `log_gamma(x)` — `ln Γ(x)` for real `x > 0`.
- `gamma(x)` — `Γ(x) = exp(log_gamma(x))`.
- `log_beta(a, b)` — `ln B(a, b)` for real `a > 0`, `b > 0`.
- `beta(a, b)` — `B(a, b) = exp(log_beta(a, b))`.

Run the tests:

```
PYTHONPATH=src python -m unittest discover -s tests
```

## Design notes

The window stores values eagerly rather than keeping running aggregates. Running
sums drift with floating point over long streams, and recomputing from a small
buffer is cheap enough that the drift is not worth the speed.

## Limitations

Values are coerced to floats, so very large integers lose precision. If you need
exact integer aggregates over a window, this is the wrong tool.

