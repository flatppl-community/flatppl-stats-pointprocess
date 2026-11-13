"""Independent count-law, finite-window, marked-event, and posterior oracles."""

import json
import math
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNTIME_POINTS = [[2.0, 3.0, -0.4], [0.7, 1.4, 0.8]]


def poisson_logpmf(n, mean):
    if mean == 0:
        return 0.0 if n == 0 else -math.inf
    return n * math.log(mean) - mean - math.lgamma(n + 1)


def runtime_oracle(x):
    """Return independent values and analytic Jacobian for runtime.flatppl."""
    signal, background, logit_efficiency = x
    efficiency = 1 / (1 + math.exp(-logit_efficiency))
    scale = 2 * efficiency
    scale_derivative = 2 * efficiency * (1 - efficiency)
    expected = scale * (signal + background)
    expected_derivative = [scale, scale, scale_derivative * (signal + background)]

    def intensities(signal_shape, background_shape):
        values, derivatives = [], []
        for fs, fb in zip(signal_shape, background_shape):
            intrinsic = signal * fs + background * fb
            values.append(scale * intrinsic)
            derivatives.append([scale * fs, scale * fb, scale_derivative * intrinsic])
        return values, derivatives

    def score(values, derivatives):
        # Count-law times conditional iid positions, then change reference measure.
        count = len(values)
        conditional = sum(math.log(value / expected) for value in values)
        labelled = poisson_logpmf(count, expected) + conditional
        configuration = labelled + math.lgamma(count + 1)
        log_derivative = [
            sum(d[j] / value for d, value in zip(derivatives, values)) for j in range(3)
        ]
        configuration_derivative = [
            log_derivative[j] - expected_derivative[j] for j in range(3)
        ]
        conditional_derivative = [
            log_derivative[j] - count * expected_derivative[j] / expected
            for j in range(3)
        ]
        return (
            configuration,
            labelled,
            conditional,
            configuration_derivative,
            conditional_derivative,
        )

    unmarked, unmarked_derivative = intensities([0.2, 0.6, 1.0], [1.0, 1.0, 1.0])
    marked, marked_derivative = intensities([0.16, 0.12, 0.8], [0.3, 0.7, 0.3])
    configuration, labelled, conditional, derivative, conditional_derivative = score(
        unmarked, unmarked_derivative
    )
    marked_score = score(marked, marked_derivative)
    return (
        [
            unmarked[0],
            expected,
            configuration,
            labelled,
            conditional,
            marked[1],
            marked_score[0],
        ],
        [
            unmarked_derivative[0],
            expected_derivative,
            derivative,
            derivative,
            conditional_derivative,
            marked_derivative[1],
            marked_score[3],
        ],
    )


def decode(value):
    if isinstance(value, str):
        return float(value)
    if isinstance(value, dict) and "data" in value:
        data = value["data"]
        return list(data.values()) if isinstance(data, dict) else data
    return value


def evaluate(request):
    result = subprocess.run(
        ["node", str(ROOT / "tests/evaluate.cjs")],
        input=json.dumps(request),
        text=True,
        capture_output=True,
        check=True,
    )
    return [decode(value) for value in json.loads(result.stdout)]


def close(actual, expected, tolerance=3e-12):
    if isinstance(expected, (list, tuple)):
        assert len(actual) == len(expected)
        for a, e in zip(actual, expected):
            close(a, e, tolerance)
    else:
        assert actual == expected or math.isclose(
            actual, expected, rel_tol=tolerance, abs_tol=tolerance
        ), (actual, expected)


def run():
    # Exact Poisson cell probabilities recover the continuous configuration density.
    rate = 1.2
    times = [0.2, 0.8, 1.1, 1.9]
    width = 0.01
    cell_means = [
        rate
        * (
            (t + width / 2) ** 2 / 2
            + t
            + width / 2
            - (t - width / 2) ** 2 / 2
            - (t - width / 2)
        )
        for t in times
    ]
    expected_count = 4 * rate
    cell_score = sum(poisson_logpmf(1, mean) for mean in cell_means)
    cell_score += poisson_logpmf(0, expected_count - sum(cell_means))
    cell_score -= len(times) * math.log(width)
    intensities = [rate * (1 + t) for t in times]
    vector = json.dumps(intensities)
    configuration, labelled, conditional = evaluate(
        [
            f"a.configuration_loglikelihood({vector},{expected_count})",
            f"a.labelled_loglikelihood({vector},{expected_count})",
            f"a.conditional_labelled_loglikelihood({vector},{expected_count})",
        ]
    )
    close(configuration, cell_score)
    close(labelled, poisson_logpmf(4, expected_count) + conditional)
    # Volume of ordered simplex is T^n/n!, while labelled cube has volume T^n.
    n, window, homogeneous_rate = 3, 2.0, 1.5
    constant = json.dumps([homogeneous_rate] * n)
    ordered, labelled = evaluate(
        [
            f"a.configuration_loglikelihood({constant},{window * homogeneous_rate})",
            f"a.labelled_loglikelihood({constant},{window * homogeneous_rate})",
        ]
    )
    close(
        math.exp(ordered) * window**n / math.factorial(n),
        math.exp(poisson_logpmf(n, window * homogeneous_rate)),
    )
    close(
        math.exp(labelled) * window**n,
        math.exp(poisson_logpmf(n, window * homogeneous_rate)),
    )
    # The empty process and structural zero rates retain their exact probabilities.
    close(
        evaluate(
            [
                "a.configuration_loglikelihood(fill(0.0,0),0.0)",
                "a.configuration_loglikelihood(fill(0.0,0),2.0)",
                "a.labelled_loglikelihood(fill(0.0,0),0.0)",
                "a.configuration_loglikelihood([0.0,1.0],1.0)",
                "a.labelled_loglikelihood([0.0],0.0)",
                "a.conditional_labelled_loglikelihood([0.0],1.0)",
                "a.conditional_labelled_loglikelihood(fill(0.0,0),1.0)",
            ]
        ),
        [0, -2, 0, -math.inf, -math.inf, -math.inf, 0],
    )
    # Count conditioning removes any common positive rate scale; order is immaterial.
    close(
        evaluate(
            [
                "a.conditional_labelled_loglikelihood([1.0,2.0],3.0)",
                "a.conditional_labelled_loglikelihood([7.0,14.0],21.0)",
                "a.configuration_loglikelihood([2.0,1.0],3.0)",
            ]
        ),
        [math.log(2 / 9), math.log(2 / 9), math.log(2) - 3],
    )
    # Binary mark sums preserve intensity before mark-dependent selection.
    close(
        evaluate(
            [
                "sum(a.marked_intensity([2.0,2.0],[0.3,0.7]))",
                "a.quadrature_integral(a.detected_intensity([1.0,2.0,3.0],[2.0,2.0,2.0],[0.5,0.5,0.5]),[0.5,1.0,0.5])",
            ]
        ),
        [2.0, 4.0],
    )
    # Explicit expected counts are used unchanged, even without normalized shapes.
    close(evaluate(["a.configuration_loglikelihood([2.0],7.0)"]), [math.log(2) - 7])
    source = (ROOT / "tests/runtime.flatppl").read_text()
    for x in RUNTIME_POINTS:
        fixed = source.replace("elementof(cartpow(reals, 3))", json.dumps(x))
        close(
            evaluate({"source": fixed, "bindings": ["outputs"]})[0],
            runtime_oracle(x)[0],
        )
    # Actual posterior density equals a normalized Gamma density plus exact evidence.
    log_evidence = 4 * math.log(4) - math.lgamma(5) + math.lgamma(6) - 6 * math.log(5)
    expressions = [f"logdensityof(p.Pi_post,record(rate={r}))" for r in [0.3, 1.2, 2.7]]
    posterior_source = (
        'flatppl_compat="0.1"\np=load_module("../examples/playground.flatppl")\n'
    )
    posterior_source += "\n".join(
        f"score{i}={expression}" for i, expression in enumerate(expressions)
    )
    posterior = evaluate(
        {"source": posterior_source, "bindings": [f"score{i}" for i in range(3)]}
    )
    close(
        posterior,
        [
            6 * math.log(5) - math.lgamma(6) + 5 * math.log(r) - 5 * r + log_evidence
            for r in [0.3, 1.2, 2.7]
        ],
    )
    # Integrate actual configuration scores on each ordered simplex, then sum counts.
    scores = evaluate(
        [f"a.configuration_loglikelihood(fill(1.5,{n}),3.0)" for n in range(40)]
    )
    mass = sum(
        math.exp(score) * 2**n / math.factorial(n) for n, score in enumerate(scores)
    )
    close(mass, 1.0)
    print(
        json.dumps(
            {
                "runtime_points": len(RUNTIME_POINTS),
                "runtime_outputs": 7,
                "posterior_points": 3,
                "normalization_mass": mass,
                "cell_density_error": abs(configuration - cell_score),
            }
        )
    )


if __name__ == "__main__":
    run()
