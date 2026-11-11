# Numerical scope and conventions

## Amplitudes

The integration drivers evaluate connected iterations of the cubic generator:

\[
 A_{22}^{\mathrm{it}}=-\frac12\langle f_2|T_1^2|i_2\rangle_{\mathrm{conn}},
 \qquad
 A_{23}^{\mathrm{it}}=-\frac{i}{6}\langle f_3|T_1^3|i_2\rangle_{\mathrm{conn}}.
\]

The full leading connected `2 -> 3` coefficient would also include
`i T3` and `-(T1 T2 + T2 T1)/2`, coherently before squaring. These are
not provided by tighter quadrature or smaller statistical errors.

The adopted compositional input cancels residual quartic sectors of `T2`
under that input's assumptions. This does not evaluate an independent
fourth-order NSF outgoing kernel or prove that the quintic contact vanishes.
`algebra/bch.py` accepts an independently supplied fourth-order coefficient,
subtracts lower-order BCH iterations and tests the reconstructed quintic.
Synthetic controls of that inverse problem are not new physical amplitudes.

## Measures and domain

Use `Ei=G=epsilon=1`, `dmu(p)=d^3p/(2|p|)`, with no added `(2pi)^-3` in
that measure. External mode factors are already part of the reduced
oscillator integrand. Incoming helicities are averaged over four choices;
outgoing helicities are summed. The identical-final-particle factors are
`1/2!` and `1/3!`. Spatial closure is resolved exactly and there is no
external energy-conserving delta.

`Emax` bounds the sum of final energies. The IR restriction applies to
outgoing momentum norms and the nonzero retained internal two-leaf
routings. The exact incoming-pair spatial zero mode is outside the
radiative domain in the center of mass.

There is no fixed excluded angular cone in these drivers. Cubic angular
kernels use a circular finite-part extension, with caps integrated to zero
radius and their complement included. The default partition fraction is
0.35; it is not 0.35 radians. Exactly collinear vertices have no assigned
ordinary point value. A partition test does not test a different
distributional extension.

## Means and physical variance

For the reference variance data, the note's moments are `P`, `W`, `Q`, `C`:
normalization, first energy moment, raw second moment and second moment
centered at the incoming energy. The older integration code uses the
name `Q` for its **first** moment, which is `W` in the variance report.

\[
 r_E=(W/P)-1,\qquad v=\mathrm{Var}(E_f/E_i)=C/P-r_E^2.
\]

The reduction computes centered moments directly, cross-checks the raw
second-moment formula, and verifies the within/between-channel decomposition.
Folded `2 -> 2` nodes combine energies `u` and `1/u`; treating every
contribution as having energy `u` gives an incorrect second moment.

## Accuracy

The direct small-cutoff runs use eight independent scrambled Sobol replicas,
256 low-energy and 128 high-energy points per replica, and angular order 16.
The retained `1e-5` reference uses four replicas and a recorded
inclusion-probability angular correction; it must not be assigned the
same sampling budget. The covariance of the moments is estimated jointly
between replicas, with independent-stratum covariances added.

Errors on a mean or a variance are estimator standard errors, not physical
energy widths. They are not generally obtained by dividing the physical
variance by the number of importance-sampling points. The deterministic
`2 -> 2` moments are held fixed in the statistical covariance; radial
quadrature changes are reported separately. Pointwise angular tests do not
provide a uniform integral-error bound.

The standalone drivers use preassigned sample counts, geometry-dependent
panels and selected order-doubling checks, not adaptive stopping at a
target statistical tolerance. Their direct acceptance tests use relative
changes below `1e-3` for selected angular weights and the radial `I22` check.

At IR `1e-8`, the paired reference gives physical variance about `4.02e-5`
at `Emax/Ei=2024` and `4.11e-4` at 18096. Mean stability alone does not
establish second-moment or ultraviolet convergence. All conclusions remain
conditional on the retained amplitudes, finite domains and angular prescription.
