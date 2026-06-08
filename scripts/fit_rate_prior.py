#!/usr/bin/env python3
"""
Fit a log-normal prior from posterior summary statistics.

Given either:
- mean + 95% interval (lower, upper), or
- median + 95% interval (lower, upper),

this script fits (mu, sigma) for a log-normal distribution by minimizing
relative squared error against the provided summary values.
"""

from __future__ import annotations

import argparse
import math
import sys
from statistics import NormalDist
from typing import Dict, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit log-normal parameters (mu, sigma) from posterior summary stats "
            "(mean or median + interval bounds)."
        )
    )

    central = parser.add_mutually_exclusive_group(required=True)
    central.add_argument(
        "--mean",
        type=float,
        help="Posterior mean in real space (must be > 0)",
    )
    central.add_argument(
        "--median",
        type=float,
        help="Posterior median in real space (must be > 0)",
    )

    parser.add_argument(
        "--lower",
        type=float,
        required=True,
        help="Lower bound of the posterior interval (must be > 0)",
    )
    parser.add_argument(
        "--upper",
        type=float,
        required=True,
        help="Upper bound of the posterior interval (must be > 0 and > lower)",
    )
    parser.add_argument(
        "--interval",
        type=float,
        default=0.95,
        help="Central interval mass, default 0.95 (e.g. 0.95 for 95%%)",
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    for name in ("mean", "median", "lower", "upper"):
        value = getattr(args, name, None)
        if value is not None and value <= 0:
            raise ValueError(f"{name} must be > 0")

    if args.lower >= args.upper:
        raise ValueError("lower must be < upper")

    if not (0.0 < args.interval < 1.0):
        raise ValueError("interval must be in (0, 1)")


def lognormal_stats(mu: float, sigma: float, z_lower: float, z_upper: float) -> Dict[str, float]:
    median = math.exp(mu)
    mean = math.exp(mu + 0.5 * sigma * sigma)
    variance = (math.exp(sigma * sigma) - 1.0) * math.exp(2.0 * mu + sigma * sigma)
    stdev = math.sqrt(variance)
    lower = math.exp(mu + sigma * z_lower)
    upper = math.exp(mu + sigma * z_upper)

    return {
        "mean": mean,
        "median": median,
        "stdev": stdev,
        "lower": lower,
        "upper": upper,
    }


def objective(
    mu: float,
    log_sigma: float,
    central_name: str,
    central_target: float,
    lower_target: float,
    upper_target: float,
    z_lower: float,
    z_upper: float,
) -> float:
    sigma = math.exp(log_sigma)
    fitted = lognormal_stats(mu, sigma, z_lower, z_upper)

    # Relative residuals keep scale balanced for rates of different magnitudes.
    r_central = (fitted[central_name] - central_target) / central_target
    r_lower = (fitted["lower"] - lower_target) / lower_target
    r_upper = (fitted["upper"] - upper_target) / upper_target

    return r_central * r_central + r_lower * r_lower + r_upper * r_upper


def fit_lognormal(
    central_name: str,
    central_target: float,
    lower_target: float,
    upper_target: float,
    interval: float,
) -> Tuple[float, float, float]:
    alpha = (1.0 - interval) / 2.0
    nd = NormalDist()
    z_lower = nd.inv_cdf(alpha)
    z_upper = nd.inv_cdf(1.0 - alpha)

    ln_lower = math.log(lower_target)
    ln_upper = math.log(upper_target)

    sigma0 = (ln_upper - ln_lower) / (z_upper - z_lower)
    sigma0 = max(sigma0, 1e-9)

    if central_name == "median":
        mu0 = math.log(central_target)
    else:
        mu0 = math.log(central_target) - 0.5 * sigma0 * sigma0

    x0 = [mu0, math.log(sigma0)]

    def f(x):
        return objective(
            x[0],
            x[1],
            central_name,
            central_target,
            lower_target,
            upper_target,
            z_lower,
            z_upper,
        )

    best = x0[:]
    best_obj = f(best)

    # Coordinate pattern search: robust and dependency-free for this smooth 2D fit.
    steps = [0.5, 0.3]
    for _ in range(300):
        improved = False

        for idx in range(2):
            for direction in (-1.0, 1.0):
                cand = best[:]
                cand[idx] += direction * steps[idx]
                cand_obj = f(cand)
                if cand_obj < best_obj:
                    best, best_obj = cand, cand_obj
                    improved = True

        # Try diagonal moves as well to avoid axis-only stagnation.
        for d0 in (-1.0, 1.0):
            for d1 in (-1.0, 1.0):
                cand = [best[0] + d0 * steps[0], best[1] + d1 * steps[1]]
                cand_obj = f(cand)
                if cand_obj < best_obj:
                    best, best_obj = cand, cand_obj
                    improved = True

        if not improved:
            steps[0] *= 0.7
            steps[1] *= 0.7
            if max(steps) < 1e-7:
                break

    mu = best[0]
    sigma = math.exp(best[1])
    return mu, sigma, best_obj


def fmt(x: float) -> str:
    return f"{x:.10g}"


def main() -> int:
    args = parse_args()

    try:
        validate_args(args)
    except ValueError as exc:
        print(f"Input error: {exc}", file=sys.stderr)
        return 2

    central_name = "mean" if args.mean is not None else "median"
    central_value = args.mean if args.mean is not None else args.median

    mu, sigma, loss = fit_lognormal(
        central_name=central_name,
        central_target=central_value,
        lower_target=args.lower,
        upper_target=args.upper,
        interval=args.interval,
    )

    alpha = (1.0 - args.interval) / 2.0
    nd = NormalDist()
    z_lower = nd.inv_cdf(alpha)
    z_upper = nd.inv_cdf(1.0 - alpha)

    fitted = lognormal_stats(mu, sigma, z_lower, z_upper)

    print("Fitted log-normal parameters")
    print(f"mu (log-space):    {fmt(mu)}")
    print(f"sigma (log-space): {fmt(sigma)}")
    print()

    print("Real-space moments of fitted distribution")
    print(f"mean:  {fmt(fitted['mean'])}")
    print(f"stdev: {fmt(fitted['stdev'])}")
    print()

    print(f"Fitted summary for central {args.interval * 100:.2f}% interval")
    print(f"mean:   {fmt(fitted['mean'])}")
    print(f"median: {fmt(fitted['median'])}")
    print(f"lower:  {fmt(fitted['lower'])}")
    print(f"upper:  {fmt(fitted['upper'])}")
    print()

    print("Deviations from provided posterior summary")
    provided = {
        central_name: central_value,
        "lower": args.lower,
        "upper": args.upper,
    }

    for key in (central_name, "lower", "upper"):
        target = provided[key]
        fit = fitted[key]
        abs_dev = fit - target
        rel_dev_pct = 100.0 * abs_dev / target
        print(
            f"{key:>6}: target={fmt(target)} fitted={fmt(fit)} "
            f"abs_dev={fmt(abs_dev)} rel_dev={rel_dev_pct:.6f}%"
        )

    print()
    print(f"Objective (sum of squared relative residuals): {loss:.12g}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
