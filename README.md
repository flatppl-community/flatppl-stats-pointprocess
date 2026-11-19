# flatppl-stats-pointprocess

[Function reference](docs/reference.md) · [Live example](https://live.flatppl.org/#model=https%3A%2F%2Fraw.githubusercontent.com%2Fflatppl-community%2Fflatppl-stats-pointprocess%2Fmain%2Fexamples%2Fplayground.flatppl&target=Pi_post)

Finite-window Poisson process scores and intensities for FlatPPL. Library version **0.1**.

For event intensities $\lambda_i$ and integrated observed intensity $\Lambda$ over the same window,

$$\ell_{\rm configuration}=\sum_{i=1}^n\log\lambda_i-\Lambda,\qquad
\ell_{\rm labelled}=\ell_{\rm configuration}-\log(n!).$$

```flatppl
flatppl_compat = "0.1"
p = load_module("main.flatppl")
intensity = p.signal_background_intensity([0.4, 0.8], [0.5, 0.5], 3.0, 2.0)
score = p.configuration_loglikelihood(intensity, 5.0)
```

Likelihood inputs are a finite nonnegative intensity vector and scalar $\Lambda\geq0$; an event at zero intensity has score $-\infty$. Conditional labelled scores subtract $n\log\Lambda$ instead of $\Lambda$ and require $\Lambda>0$. Intensity helpers combine matching vectors or scalars; quadrature takes matching node-value and weight vectors.

The [rate example](examples/playground.flatppl) has an analytic Gamma posterior. See the reference for factorial conventions, detection/mark normalization, generation and validation.
