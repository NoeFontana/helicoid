"""Function ids. A new id is an `evaluate` in a math module plus one line in `FUNCTIONS`."""

from collections.abc import Callable
from dataclasses import dataclass

from . import check, coeff
from .strata import R_STRATA, SCALAR_THETA_STRATA, Stratum


@dataclass(frozen=True)
class FunctionSpec:
    name: str
    strata: tuple[Stratum, ...]
    inputs: Callable[[Stratum], list[dict]]  # exact binary64 inputs: float or list[float] values
    evaluate: Callable[[dict], dict]  # definition at the current mp.dps: mpf or list[mpf] values
    check: Callable[[dict, dict], None] | None = None  # (inputs, out), raises on a disagreement


FUNCTIONS: dict[str, FunctionSpec] = {
    spec.name: spec
    for spec in (
        *(
            FunctionSpec(
                f"coeff_{c}",
                SCALAR_THETA_STRATA,
                coeff.theta_inputs,
                coeff.evaluator(c),
                check.coefficient(c),
            )
            for c in "kabcde"
        ),
        FunctionSpec("coeff_r", R_STRATA, coeff.r_inputs, coeff.r, check.coefficient("r")),
    )
}
