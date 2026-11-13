# flatppl-stats-pointprocess

Finite-window Poisson process likelihoods, signal/background intensities, detection exposure, and marked events for FlatPPL 0.1. Library version: **0.1**.

```flatppl
flatppl_compat = "0.1"
p = load_module("main.flatppl")
intensity = p.signal_background_intensity([0.4, 0.8], [0.5, 0.5], 3.0, 2.0)
score = p.configuration_loglikelihood(intensity, 5.0)
```

The caller supplies the integrated **observed** intensity over the finite observation window. An expected count is used as given: these helpers never integrate or normalize it. The [reference](docs/reference.md) distinguishes unordered configurations, ordered continuous times, and labelled points, including their factorial factors.

The [posterior example](examples/playground.flatppl) infers one event-rate parameter with an analytic Gamma posterior. The known event shape makes the count sufficient for rate, so the example uses its exact Poisson likelihood. It also exposes the full unbinned log likelihood. The library does not construct a variable-length sampling law. FlatPPL's built-in `PoissonProcess` provides that law where an engine supports it; see the reference for the distinction and a construction example.

Run the independent mathematical checks with Python's standard library and a current Node.js with TypeScript support:

```sh
FLATPPL_JS_DIR=/path/to/flatppl-js python tests/check_math.py
```

The checks use the real JS module loader, independent Poisson count factorizations, finite-window integrals, and an analytic posterior. [The runtime fixture](tests/runtime.flatppl) and `runtime_oracle`/`RUNTIME_POINTS` in [the checker](tests/check_math.py) provide seven outputs for compiled value and derivative checks.
