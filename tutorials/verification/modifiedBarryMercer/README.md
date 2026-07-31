# Modified Barry--Mercer checkerboard benchmark

This case reproduces the undrained Barry--Mercer modification from Section
5.3 of Aronson et al., *Pressure-stabilized fixed-stress iterative solutions
of compositional poromechanics*, CMAME 427 (2024), 117008.

The paper uses a unit square, a 10-by-10 Cartesian grid, `E = 10 kPa`,
`nu = 0.2`, an incompressible unit-viscosity fluid, isotropic intrinsic
permeability `1e-12 m2`, and ten time steps of 10 s.  A source
`sin(pi*t/100)` is applied to the cell centred at `(0.35 0.15)`.

Run all three variants with:

```sh
./Allrun
```

The script creates `runs/unstabilized`, `runs/explicit`, and `runs/implicit`,
runs them, creates a `.foam` marker in each case for ParaView, and writes
`runs/checkerboardMetrics.csv`. The metric is the norm of the normalized
projection of the final pressure field onto the three nonconstant Cartesian
parity modes `(-1)^i`, `(-1)^j`, and `(-1)^(i+j)`. The two stabilized cases
should reduce this metric, and their pressure fields should remain close to
one another.

The fixed-stress scales are selected independently from a first-step sweep
followed by a ten-step neighbor check. The runner uses `5.0` for the explicit
stabilization and `5.5` for the implicit stabilization. The unstabilized
negative control retains `4.0`; it reaches the 1000-iteration cap instead of
converging and advances with a checkerboard-dominated pressure field.

The first-step sweep found that explicit scales `3.5` and `3.75`, and implicit
scales from `2.0` through `3.0`, all reached the iteration cap. The local
first-step minima were `5.1` for explicit and `5.5` for implicit. Over the full
ten-step history, explicit `5.0` was slightly better than `5.1` and `5.2`, with
1524 total outer iterations. Implicit `5.5` remained the local optimum, with
1830 total outer iterations. Unstabilized scales from `2.0` through `6.0` all
reached the cap and therefore do not represent converged solutions.

This benchmark is deliberately nearly undrained.  It tests the spatial
checkerboard response, not agreement with the drained analytical
Barry--Mercer solution.  Stabilization is disabled by selecting `poroSolid`;
the other two variants select `poroSolidStab` and differ only in
`stabilizationType`.

## Unstabilized tau--Fo sweep

Run the unstabilized parameter sweep with:

```sh
./AllrunTauFo
```

Here `tau = M_Biot/K_vol` is the coupling strength and `Fo` is the
consolidation Fourier number. The runner writes every case, including a
ParaView `.foam` marker, below `tauFoRuns/`. Its summary files are
`tauFoSweep.csv` and `tauFoSurface.png`.

The script retains two distinct contraction diagnostics. The full coupled
factor `rho` comes from adjacent RMS envelopes of the dimensionless update
`max(maxDeltaP/pTolerance, maxDeltaD/DTolerance)`. The checkerboard factor
`rho_cb` comes from the iteration-to-iteration increment of the pressure
projection onto the three-dimensional Cartesian parity subspace spanned by
`(-1)^i`, `(-1)^j`, and `(-1)^(i+j)`. These three modes span the nonconstant
part of the four parity classes decoupled by the wide Cartesian stencil.
Projected increments are used because the source can force a nonzero final
subspace amplitude: for an affine iteration they obey
`delta(a)_(k+1) = A_cb delta(a)_k`. RMS blocks of their Euclidean norms
measure the dominant subspace contraction without being masked by sign
alternation. The fit stops when the stabilized modes reach their numerical
noise floor.

The surface has three panels: final checkerboard-subspace/source amplitude,
full coupled `rho`, and checkerboard-subspace `rho_cb`. The rate panels use raw
values on a `0.7`--`1.3` diverging scale centred at one. No interpolated
contours are drawn: the sampled values and the uncertainty-aware matched-scale
plots are used to classify growth. A localized source can have a nonzero final
projection onto the parity subspace even when `rho_cb` is well below one.

The CSV retains `contractionFactor`, `checkerboardContractionFactor`,
`envelopeContracting`, `monotone`, and `converged` so
that marginal cases can be diagnosed instead of silently discarded. The
last of these additionally requires the pressure and displacement update
tolerances to be met before the cap.

Run the corresponding implicit compact-stencil stabilization sweep with:

```sh
./AllrunTauFoImplicit
```

It uses the same grid and diagnostics, selects `poroSolidStab` with
`stabilizationType implicit`, and writes to `tauFoRunsImplicit/`. The default
implicit fixed-stress scale is `5.5` and the default compact-stencil
stabilization factor is `1.0`; these can be overridden with
`--fixed-stress-scale` and `--stab-factor`.

Run the primary matched scale-`2.5` comparison with:

```sh
./AllrunTauFoUnstabilizedScale2p5
./AllrunTauFoImplicitScale2p5
```

After both sweeps exist, regenerate their surfaces and create the joint
rate-regime, final-strength, convergence, iteration-history, and pressure-field
plots without rerunning the solver:

```sh
python3 plot_matched_scale.py
```

The derived plots and CSV are written below `matchedScale2p5Plots/`. The rate
classification uses linear `rho` and `rho_cb` axes. Its parameter-space regime
map uses multi-window pressure-update, displacement-update, and checkerboard
rates to distinguish robust contraction, robust growth, and unresolved
near-neutral behavior. No interpolated stability contours are drawn. The
original sweep surfaces show the raw contraction rates on `0.7`--`1.3`
diverging colour scales centred at one rather than their logarithms.

To extend only capped points that were not already robustly checkerboard
growing, run:

```sh
python3 extend_matched_scale.py --max-iterations 1200 --jobs 4
python3 plot_matched_scale.py
```

The extension writes below `matchedScale2p5Extended/`; the plotting script
automatically prefers those longer histories where present. The saturated
fluid diagnostic records the pre- and post-solve pressure-equation defects and
the pressure linear-solver residual when `pressureEquationDiagnostic true` is
set in `poroFluidCoeffs`.

The outputs are isolated in `tauFoRunsUnstabilizedScale2p5/` and
`tauFoRunsImplicitScale2p5/`. Since the fixed-stress coefficient is identical,
the comparison shows the effect of the compact checkerboard stabilization
without conflating it with fixed-stress over-stabilization.

With the three-mode diagnostic, both variants reach the absolute tolerances
at 45 of 81 points. Extending 21 unstabilized and 36 implicit capped points
from 300 to 1200 outer iterations does not add any converged points: their
late RMS envelopes remain bounded and nearly periodic instead of trending
toward the absolute tolerances. The pressure equations themselves remain
accurately solved, with a maximum post-solve RMS defect of about `2.1e-13`.

The uncertainty-aware classification leaves the unstabilized grid with 41
robustly contracting points, 15 robustly both-growing points, and 25 unresolved
near-neutral points. At all 15 growing points the full and checkerboard rates
grow together; there is no robust checkerboard-only growth region. The
implicit formulation has 42 robustly contracting and 39 unresolved points,
with no robust growth. Its checkerboard subspace contracts at 78 points and is
unresolved at three; none grows. Stabilization limits the maximum
checkerboard/source amplitude to 1.71%, compared with 146 in the unstabilized
grid, but at fixed-stress scale `2.5` it does not remove the separate
coupling-level periodic plateau.

On the default 9-by-9 grid, scale `2.5` reaches the absolute tolerances at
45 of 81 points and is robustly monotone at 33. Scale `5.5` reaches the
absolute tolerances at all 81 points, has `rho < 1` at all 81 points, and is
robustly monotone at 39. All scale-`5.5` checkerboard amplitudes remain below
1% of the source-induced pressure scale. Its nearly undrained, low-`Fo`
corner still contracts non-monotonically; this is retained as a CSV
diagnostic but does not redefine the uncertainty-aware classification.

For a controlled convergence comparison at the same fixed-stress scale, run:

```sh
./AllrunTauFoUnstabilizedScale5p5
```

This keeps stabilization disabled but uses scale `5.5`, writing exclusively
to `tauFoRunsUnstabilizedScale5p5/`. It can therefore be compared directly
with `tauFoRunsImplicit/` without deleting or changing the original
unstabilized scale-`2.0` sweep in `tauFoRuns/`.

In this matched-scale comparison, the unstabilized formulation reaches the
absolute tolerances at 51 of 81 points, whereas the implicit formulation
reaches them at all 81. Only 31 unstabilized points have a checkerboard below
1% of the source scale, compared with all 81 implicit points. At
`tau = 1e6` and `Fo = 1e-5`, the unstabilized iteration is monotone but nearly
neutral (`rho = 0.99977`), remains unconverged after 300 iterations, and has a
checkerboard/source ratio of `1.09`. The corresponding implicit point has
`rho = 0.933`, converges in 179 iterations, and has a ratio of `0.00288`.
