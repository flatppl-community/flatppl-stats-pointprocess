# Point processes reference

Library version **0.1** · [Package overview](../README.md)

Community library version **0.1**, loaded with `load_module("main.flatppl")`. All functions return numeric values. No prior, event window, integration algorithm, normalization, or event selection is implicit.

| Function | Arguments | Result |
|---|---|---|
| [`configuration_loglikelihood`](#configuration_loglikelihood) | `event_intensities, expected_count` | unordered configuration score; also ordered continuous-time score |
| [`labelled_loglikelihood`](#labelled_loglikelihood) | `event_intensities, expected_count` | joint count and labelled-point score |
| [`conditional_labelled_loglikelihood`](#conditional_labelled_loglikelihood) | `event_intensities, expected_count` | labelled positions conditional on their count |
| [`signal_background_intensity`](#signal_background_intensity) | `signal_density, background_density, expected_signal, expected_background` | sum of two component intensities |
| [`detected_intensity`](#detected_intensity) | `intrinsic_intensity, exposure, efficiency` | observed rate after exposure and thinning |
| [`marked_intensity`](#marked_intensity) | `spatial_intensity, mark_density` | joint point-and-mark intensity |
| [`quadrature_integral`](#quadrature_integral) | `intensity, weights` | caller-selected weighted integral |

## Event space and reference measures

Let $W$ be a fixed finite observation window, $\rho$ a fixed parameter-independent event reference measure, and $\lambda=d\mu/d\rho$ a nonnegative integrable observed intensity. The expected detected count is

```math
\Lambda=\mu(W)=\int_W\lambda(t)\,d\rho(t)<\infty.
```

Use Lebesgue measure for continuous coordinates, counting measure for discrete marks, and their product for marked continuous events. Every supplied intensity value and integral must refer to the **same observed process and window**. Window membership, sorting, array lengths, nonnegativity, and consistency with the integral are caller obligations; helpers accept evaluated intensities rather than event coordinates. They do not check these obligations or repair invalid inputs.

The configuration reference measure $\nu$ is defined for symmetric integrands by

```math
\int F\,d\nu=\sum_{n=0}^{\infty}\frac{1}{n!}\int_{W^n}F(\{t_1,\ldots,t_n\})\,d\rho^n.
```

It includes multiplicities for atomic event spaces. The Poisson configuration density is $e^{-\Lambda}\prod_i\lambda(t_i)$ and integrates to $e^{-\Lambda}\sum_n\Lambda^n/n!=1$. Event array order has no meaning under this convention. This matches FlatPPL's [Poisson processes](https://github.com/flatppl/flatppl-design/blob/main/docs/06-measure-algebra.md#poisson-processes).

For **strictly ordered continuous event times** $t_1\lt\cdots\lt t_n$, the reference measure is counting measure in $n$ times ordinary Lebesgue measure restricted to the ordered simplex. The same density applies: the simplex contributes the $1/n!$ volume. Sorting iid continuous points creates this representation. This statement requires non-atomic time coordinates and is not a formula for tied discrete times.

For **a random count and iid labelled points**, the reference measure is counting in $n$ times $\rho^n$ on all labelled tuples. Its density has an additional $1/n!$. A labelled array and an unordered configuration therefore cannot share the same absolute density convention. At fixed observed $n$ this factor is parameter-independent; retaining it matters for normalization, evidence, and comparing different counts or observation conventions.

### `configuration_loglikelihood`

**`configuration_loglikelihood(event_intensities, expected_count)`** returns

$$\ell_{\rm config}=-\Lambda+\sum_{i=1}^{n}\log\lambda_i.$$

`event_intensities` is a vector of length $n$ containing finite nonnegative rates evaluated at the observed events, and `expected_count` is finite and nonnegative. An empty vector gives $-\Lambda$, including $0$ when $\Lambda=0$. Any observed event of zero intensity gives $-\infty$. Negative rates are outside the contract. The function neither divides by $n!$ nor by $\Lambda^n$.

### `labelled_loglikelihood`

**`labelled_loglikelihood(event_intensities, expected_count)`** returns

$$\ell_{\rm labelled}=\ell_{\rm config}-\log\Gamma(n+1).$$

Domains are those of `configuration_loglikelihood`. Derivation: multiply the Poisson count mass $e^{-\Lambda}\Lambda^n/n!$ by $n$ iid labelled position densities $\prod_i\lambda_i/\Lambda$. The limiting empty zero-intensity process has mass one.

### `conditional_labelled_loglikelihood`

**`conditional_labelled_loglikelihood(event_intensities, expected_count)`** returns

$$\ell_{\rm positions\mid n}=\sum_i\log\lambda_i-n\log\Lambda.$$

Requires finite $\Lambda>0$ and the same vector domain as above. This is a density with respect to $\rho^n$, **conditional on the observed count**. It discards all information carried by the total count and is invariant under a common positive rescaling of $\lambda$ and $\Lambda$. For ordered continuous positions conditional on $n$, add $\log\Gamma(n+1)$. Conditioning a zero-intensity process on a positive count is undefined.

### `signal_background_intensity`

**`signal_background_intensity(signal_density, background_density, expected_signal, expected_background)`** returns the pointwise intensity

$$\lambda(t)=S f_s(t)+B f_b(t).$$

The shapes are matching vectors of nonnegative values with the same reference measure. $S,B\geq0$ are finite scalar coefficients. If each shape integrates to one **over the observation window**, $S$ and $B$ are expected counts and $\Lambda=S+B$. Otherwise the integral is $S\int_W f_s+B\int_W f_b$; the function does not normalize the shapes. A full-line probability density truncated to a finite window generally needs a separate explicit normalization if its coefficient is meant to denote a window yield.

### `detected_intensity`

**`detected_intensity(intrinsic_intensity, exposure, efficiency)`** returns

$$\lambda_{\rm obs}(t)=r(t)E(t)\epsilon(t).$$

Inputs may be scalars or matching vectors, with finite $r,E\geq0$ and $0\leq\epsilon\leq1$; scalar inputs broadcast. Exposure carries the units converting the intrinsic rate into an event intensity, while efficiency is a detection probability. Integrate this **observed** intensity to obtain $\Lambda_{\rm obs}$. Applying a scalar efficiency to event terms while leaving $\Lambda$ unchanged is incorrect. If coefficients already denote detected yields, do not apply the same efficiency again.

### `marked_intensity`

**`marked_intensity(spatial_intensity, mark_density)`** returns

$$\lambda(t,m)=\lambda(t)q(m\mid t).$$

Inputs are finite nonnegative scalars or matching vectors evaluated at the same events. $q(\cdot\mid t)$ must integrate/sum to one in the mark reference measure; then the integrated count is unchanged. For component-specific mark laws, form each marked component first, then call `signal_background_intensity`. Multiplying the mixture by one common mark law is appropriate only when both components share that law. If selection depends on the mark, apply it before computing the observed integral; the count will generally change.

### `quadrature_integral`

**`quadrature_integral(intensity, weights)`** returns

$$\widehat\Lambda=\sum_j w_j\lambda(u_j).$$

Inputs are matching vectors. Intensities and weights are finite and nonnegative, and weights include integration units. Nodes, weights, window, accuracy, and any parameter dependence are entirely caller-controlled. This operation is exact for a supplied finite weighted sum; its approximation to a continuous integral is only as good as the supplied quadrature. It does not integrate the observed event array unless that array independently constitutes a valid quadrature rule.

## Posterior and generative use

[The rate example](../examples/playground.flatppl) uses $\lambda(t)=r(1+t)$ on $[0,2]$, with observed times $[0.2,0.8,1.1,1.9]$, so $\Lambda=4r$. With the shape/rate prior $r\sim\mathrm{Gamma}(2,1)$,

$$\pi(r\mid t)\propto r^{2-1+4}e^{-(1+4)r},\qquad r\mid t\sim\mathrm{Gamma}(6,5).$$

The mean is $1.2$ and variance $0.24$. Since the shape is known, the observed count is sufficient for $r$. The example therefore uses $N\sim\mathrm{Poisson}(4r)$ and the ordinary `bayesupdate` pipeline. The trapezoidal integral of the linear shape is exact. The full configuration log likelihood differs from the count log likelihood by $\sum_i\log(1+t_i)+\log(4!)-4\log4$, independent of $r$. The example exposes it as `unbinned_loglikelihood` too. For a shape parameter inferred from event positions, retain the position terms, for example by weighting a proper prior with `logweighted`; the JS viewer currently uses IS for that raw weighted measure.

For genuine variable-length generation, use FlatPPL's built-in law where supported:

```flatppl
flatppl_compat = "0.1"
events ~ PoissonProcess(superpose(
    weighted(3.0, Uniform(interval(0, 2))),
    weighted(2.0, Uniform(interval(0, 2)))
))
```

This is a random count followed by exactly that many points, represented by the built-in ragged process law. It is not an independent random count paired with a separately fixed-length array. Current JS supports scalar finite-window component shapes; current StableHLO does not lower `PoissonProcess`. The numeric helpers avoid that operator and provide the compiled inference path. This library provides no replacement variable-length constructor or simulator.

## Validation

Run the independent mathematical checks with Python's standard library and a current Node.js with TypeScript support:

```sh
FLATPPL_JS_DIR=/path/to/flatppl-js python tests/check_math.py
```

The checks use the real JS module loader, independent Poisson count factorizations, finite-window integrals, and an analytic posterior. [The runtime fixture](../tests/runtime.flatppl) and `runtime_oracle`/`RUNTIME_POINTS` in [the checker](../tests/check_math.py) provide seven outputs for compiled value and derivative checks.
