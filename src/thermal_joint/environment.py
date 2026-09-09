"""Thermal environment: reference state plus cold and hot temperature extremes."""

from __future__ import annotations

from dataclasses import dataclass

from ._validation import require_finite

__all__ = ["ThermalEnvironment"]


@dataclass(frozen=True)
class ThermalEnvironment:
    """Reference, cold and hot temperatures of a thermal excursion.

    Attributes
    ----------
    reference_temperature:
        Stress-free reference temperature ``T_ref``.
    cold_temperature:
        Cold extreme ``T_cold``.
    hot_temperature:
        Hot extreme ``T_hot``.

    All three must be finite and given on the same scale (all K or all degC),
    and must satisfy ``cold <= reference <= hot``. Equality is allowed, which
    makes the corresponding excursion a zero-``dT`` case.

    Only *differences* enter the mechanics, so K and degC give identical
    results as long as one scale is used consistently.
    """

    reference_temperature: float
    cold_temperature: float
    hot_temperature: float

    def __post_init__(self) -> None:
        t_ref = require_finite(self.reference_temperature, "reference_temperature")
        t_cold = require_finite(self.cold_temperature, "cold_temperature")
        t_hot = require_finite(self.hot_temperature, "hot_temperature")
        if not (t_cold <= t_ref <= t_hot):
            raise ValueError(
                "temperatures must satisfy cold <= reference <= hot, got "
                f"cold={t_cold!r}, reference={t_ref!r}, hot={t_hot!r}"
            )
        object.__setattr__(self, "reference_temperature", t_ref)
        object.__setattr__(self, "cold_temperature", t_cold)
        object.__setattr__(self, "hot_temperature", t_hot)

    @property
    def delta_temperature_cold(self) -> float:
        """``T_cold - T_ref``; zero or negative."""
        return self.cold_temperature - self.reference_temperature

    @property
    def delta_temperature_hot(self) -> float:
        """``T_hot - T_ref``; zero or positive."""
        return self.hot_temperature - self.reference_temperature
