"""Function ids. A new id is an `evaluate` in a math module plus one line in `FUNCTIONS`."""

from collections.abc import Callable
from dataclasses import dataclass

from . import check, check_se2, check_sen3, coeff, se2, sen3, so2, so3
from .strata import (
    QUAT_STRATA,
    R_STRATA,
    SCALAR_THETA_STRATA,
    SE2_STRATA,
    SEN3_QUAT_STRATA,
    SEN3_STRATA,
    Stratum,
)


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
        FunctionSpec("so2_exp", SCALAR_THETA_STRATA, so2.theta_inputs, so2.exp, check_se2.so2_exp),
        FunctionSpec("so2_log", SCALAR_THETA_STRATA, so2.z_inputs, so2.log, check_se2.so2_log),
        FunctionSpec("se2_exp", SE2_STRATA, se2.tau_inputs, se2.exp, check_se2.se2_exp),
        FunctionSpec("se2_log", SE2_STRATA, se2.x_inputs, se2.log, check_se2.se2_log),
        FunctionSpec("se2_ad", SE2_STRATA, se2.x_inputs, se2.ad_of, check_se2.se2_ad),
        *(
            FunctionSpec(
                f"se2_{name}",
                SE2_STRATA,
                se2.tau_inputs,
                se2.jacobian(name),
                check_se2.jacobian(name),
            )
            for name in ("jr", "jl", "jr_inv", "jl_inv")
        ),
        *(
            spec
            for n in (1, 2, 3)
            for spec in (
                FunctionSpec(
                    f"sen3_exp_n{n}",
                    SEN3_STRATA,
                    sen3.tau_inputs(n),
                    sen3.exp(n),
                    check_sen3.exp(n),
                ),
                FunctionSpec(
                    f"sen3_log_n{n}",
                    SEN3_QUAT_STRATA,
                    sen3.x_inputs(n, both_signs=True),
                    sen3.log(n),
                    check_sen3.log(n),
                ),
                FunctionSpec(
                    f"sen3_ad_n{n}",
                    SEN3_QUAT_STRATA,
                    sen3.x_inputs(n, both_signs=False),
                    sen3.ad_of(n),
                    check_sen3.ad(n),
                ),
                *(
                    FunctionSpec(
                        f"sen3_{name}_n{n}",
                        SEN3_STRATA,
                        sen3.tau_inputs(n),
                        sen3.jacobian(n, name),
                        check_sen3.jacobian(n, name),
                    )
                    for name in ("jr", "jl", "jr_inv", "jl_inv")
                ),
            )
        ),
    )
}
