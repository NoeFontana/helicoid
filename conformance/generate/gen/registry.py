"""Function ids. A new id is an `evaluate` in a math module plus one line in `FUNCTIONS`."""

from collections.abc import Callable
from dataclasses import dataclass

from . import check, coeff, so3
from .strata import QUAT_STRATA, R_STRATA, SCALAR_THETA_STRATA, Stratum


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
        FunctionSpec("so3_exp", SCALAR_THETA_STRATA, so3.phi_inputs, so3.exp, check.so3_exp),
        FunctionSpec("so3_log", QUAT_STRATA, so3.log_inputs, so3.log, check.so3_log),
        FunctionSpec("so3_act", QUAT_STRATA, so3.act_inputs, so3.act, check.so3_act),
        FunctionSpec(
            "so3_from_matrix",
            QUAT_STRATA,
            so3.matrix_inputs,
            so3.from_matrix,
            check.so3_from_matrix,
        ),
        *(
            FunctionSpec(
                f"so3_{name}",
                SCALAR_THETA_STRATA,
                so3.phi_inputs,
                so3.jacobian(name),
                check.so3_jacobian(name),
            )
            for name in ("jr", "jl", "jr_inv", "jl_inv")
        ),
    )
}
