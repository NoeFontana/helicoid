"""Function ids. A new id is an `evaluate` in a math module plus one line in `FUNCTIONS`."""

from collections.abc import Callable
from dataclasses import dataclass

from . import coeff
from .strata import SCALAR_THETA_STRATA, Stratum


@dataclass(frozen=True)
class FunctionSpec:
    name: str
    strata: tuple[Stratum, ...]
    inputs: Callable[[Stratum], list[dict]]  # exact binary64 inputs: float or list[float] values
    evaluate: Callable[[dict], dict]  # definition at the current mp.dps: mpf or list[mpf] values


FUNCTIONS: dict[str, FunctionSpec] = {
    spec.name: spec
    for spec in (FunctionSpec("coeff_k", SCALAR_THETA_STRATA, coeff.theta_inputs, coeff.k),)
}
