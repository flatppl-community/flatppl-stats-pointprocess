# flatppl-stats-pointprocess

[Function reference](docs/reference.md) · [Live example](https://live.flatppl.org/#model=https%3A%2F%2Fraw.githubusercontent.com%2Fflatppl-community%2Fflatppl-stats-pointprocess%2Fmain%2Fexamples%2Fplayground.flatppl&target=Pi_post)

Finite-window Poisson process likelihoods, signal/background intensities, detection exposure, and marked events for FlatPPL. Library version **0.1**.

```flatppl
flatppl_compat = "0.1"
p = load_module("main.flatppl")
intensity = p.signal_background_intensity([0.4, 0.8], [0.5, 0.5], 3.0, 2.0)
score = p.configuration_loglikelihood(intensity, 5.0)
```

The caller supplies the integrated observed intensity. Helpers use it without integration or normalization. Configuration, ordered-time, and labelled-point likelihoods have distinct factorial conventions.

The [rate example](examples/playground.flatppl) has an analytic Gamma posterior. It exposes count and unbinned likelihoods. The library does not construct a variable-length sampling law.

See the reference for definitions, domains, numerical limits, and validation commands.
