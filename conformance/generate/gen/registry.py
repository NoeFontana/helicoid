"""Function ids. A new id is an `evaluate` in a math module plus one line in `FUNCTIONS`."""

from collections.abc import Callable
from dataclasses import dataclass

from . import (
    charts,
    check,
    check_charts,
    check_geodesic,
    check_linalg,
    check_real,
    check_se2,
    check_sen3,
    chol,
    coeff,
    cubic,
    eig,
    geodesic,
    mat2,
    quat,
    real,
    se2,
    sen3,
    so2,
    so3,
)
from .strata import (
    COEFF_R_STRATA,
    COEFF_STRATA,
    GEO_STRATA,
    QUAT_STRATA,
    SCALAR_THETA_STRATA,
    SE2_STRATA,
    SEN3_QUAT_STRATA,
    SEN3_STRATA,
    Drawn,
    Stratum,
    drawn_inputs,
)


@dataclass(frozen=True)
class FunctionSpec:
    name: str
    strata: tuple[Stratum | Drawn, ...]
    inputs: Callable[[Stratum], list[dict]]  # exact binary64 inputs: float or list[float] values
    evaluate: Callable[[dict], dict]  # definition at the current mp.dps: mpf or list[mpf] values
    check: Callable[[dict, dict], None] | None = None  # (inputs, out), raises on a disagreement


FUNCTIONS: dict[str, FunctionSpec] = {
    spec.name: spec
    for spec in (
        *(
            FunctionSpec(
                f"coeff_{c}",
                COEFF_STRATA,
                coeff.theta_inputs,
                coeff.evaluator(c),
                check.coefficient(c),
            )
            for c in ("k", "a", "b", "c", "d", "e", "cos_half", "alpha")
        ),
        FunctionSpec("coeff_r", COEFF_R_STRATA, coeff.r_inputs, coeff.r, check.coefficient("r")),
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
        # 0064: the strata of `so3_jl` and their `@f32` twins (0016 item 1), the one vector id
        # scored at both precisions.
        FunctionSpec("so3_gamma2", COEFF_STRATA, so3.phi_inputs, so3.gamma2, check.so3_gamma2),
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
        FunctionSpec(
            "so3_geodesic",
            GEO_STRATA,
            geodesic.so3_inputs,
            geodesic.so3_evaluate,
            check_geodesic.so3_geodesic,
        ),
        FunctionSpec(
            "se3_geodesic",
            GEO_STRATA,
            geodesic.se3_inputs,
            geodesic.se3_evaluate,
            check_geodesic.se3_geodesic,
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
        # docs/decisions/0060 decision 8: the three SE(3) charts, binary64, over the SE(3) strata.
        *(
            spec
            for chart in charts.CHARTS
            for spec in (
                FunctionSpec(
                    f"se3_{chart}_retract",
                    SEN3_STRATA,
                    charts.retract_inputs,
                    charts.retract(chart),
                    check_charts.retract(chart),
                ),
                FunctionSpec(
                    f"se3_{chart}_local",
                    SEN3_STRATA,
                    charts.local_inputs(chart),
                    charts.local(chart),
                    check_charts.local(chart),
                ),
            )
        ),
        # docs/decisions/0056: the routines D7 did not reach. Binary64 strata, then `@f32` ones.
        FunctionSpec(
            "solve_cubic", cubic.STRATA, drawn_inputs, cubic.solve_cubic, check_linalg.solve_cubic
        ),
        FunctionSpec("eig3", eig.STRATA, drawn_inputs, eig.eig3, check_linalg.eig3),
        *(
            FunctionSpec(f"chol_n{n}", chol.strata(n, False), drawn_inputs, chol.chol, check_linalg.chol)
            for n in (3, 6)
        ),
        *(
            FunctionSpec(
                f"chol_solve_n{n}",
                chol.strata(n, True),
                drawn_inputs,
                chol.chol_solve,
                check_linalg.chol_solve,
            )
            for n in (3, 6)
        ),
        # docs/decisions/0061.
        FunctionSpec(
            "mat2_inverse_adj",
            mat2.STRATA,
            drawn_inputs,
            mat2.inverse_adj,
            check_linalg.mat2_inverse_adj,
        ),
        FunctionSpec(
            "quat_renormalize",
            quat.STRATA,
            drawn_inputs,
            quat.renormalize,
            check_linalg.quat_renormalize,
        ),
        FunctionSpec("real_sqrt", real.SQRT_STRATA, drawn_inputs, real.sqrt, check_real.real_sqrt),
        FunctionSpec("real_cbrt", real.CBRT_STRATA, drawn_inputs, real.cbrt, check_real.real_cbrt),
        FunctionSpec(
            "real_sin_cos", real.SIN_COS_STRATA, drawn_inputs, real.sin_cos, check_real.real_sin_cos
        ),
        FunctionSpec("real_acos", real.ACOS_STRATA, drawn_inputs, real.acos, check_real.real_acos),
        FunctionSpec(
            "real_atan2", real.ATAN2_STRATA, drawn_inputs, real.atan2, check_real.real_atan2
        ),
        FunctionSpec("real_div", real.DIV_STRATA, drawn_inputs, real.div, check_real.real_div),
    )
}
