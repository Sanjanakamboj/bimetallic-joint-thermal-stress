"""Shared numeric input validation helpers.

All quantities in this package are SI:

    modulus  [Pa]
    CTE      [1/K]
    area     [m^2]
    stress   [Pa]
    force    [N]
    strain   [-]
"""

from __future__ import annotations

import math

__all__ = ["require_finite", "require_positive", "require_non_empty_name"]


def require_finite(value: float, name: str) -> float:
    """Return ``value`` as a float, rejecting NaN and +/- infinity."""
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:  # pragma: no cover - defensive
        raise TypeError(f"{name} must be a real number, got {value!r}") from exc
    if not math.isfinite(numeric):
        raise ValueError(f"{name} must be finite, got {numeric!r}")
    return numeric


def require_positive(value: float, name: str) -> float:
    """Return ``value`` as a finite float, rejecting zero and negatives."""
    numeric = require_finite(value, name)
    if numeric <= 0.0:
        raise ValueError(f"{name} must be strictly positive, got {numeric!r}")
    return numeric


def require_non_empty_name(value: str, name: str = "name") -> str:
    """Return a stripped, non-empty string."""
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a string, got {type(value).__name__}")
    stripped = value.strip()
    if not stripped:
        raise ValueError(f"{name} must be a non-empty string")
    return stripped
