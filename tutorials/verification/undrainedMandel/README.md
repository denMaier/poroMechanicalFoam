# Fully undrained Mandel first-increment test

> **Result provenance:** sweep directories generated with the former static
> top-displacement value are invalid because that setup produced a zero
> old-to-current displacement increment on the top patch.  Do not reuse those
> directories.  The runner now writes a `displacementSeries` ramp for every
> case; a runtime check confirmed that it retains the zero old boundary value
> and gives the prescribed nonzero increment at the first time step.

This test starts from `D=0` and `p_rgh=0`, applies a smooth frictionless
prescribed displacement during the first time step, and makes every
hydraulic boundary impermeable.  It is therefore a test of the discrete
undrained pressure response without the drained-edge boundary layer of the
classical Mandel problem.

The top displacement is ramped from zero at `t=0` to
`ubar*(1 + 0.1*cos(pi*x/a))` at the case-specific `dt`, then held constant.
The time series retains the correct old-time boundary value, while the
optional per-face `spatialScale` on `fixedDisplacementZeroShear` supplies only
a resolved, domain-scale asymmetry.  Use `--load-variation 0` for the uniform
control.

The thesis coupling strength is used directly:

    tau = M/K_vol,
    Ss  = 1/M = 1/(tau K_vol).

The incompressible-storage limit is therefore `tau -> infinity` and
`1/M = Ss -> 0`.  The second sweep coordinate is the mesh Fourier number

    Fo = Kp dt / ((Ss + 1/K_vol) h^2),

with `Kp=k/|rho g|` and the uniform cell size `h`.  Conductivity is fixed and
the runner computes `dt` from each requested `Fo`, so `Fo -> 0` represents
the first-increment limit `dt -> 0`.

Run a small matched comparison with, for example:

    ./undrained_mandel_sweep.py --tau-values "1 3 10 30 100" \
        --fo-values "1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1" \
        --formulation unstabilized --output-directory tauFo_unstabilized --clean
    ./undrained_mandel_sweep.py --tau-values "1 3 10 30 100" \
        --fo-values "1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1" \
        --formulation implicit --stab-factor 0.5 \
        --output-directory tauFo_stab0p5 --clean
    ./undrained_mandel_sweep.py --tau-values "1 3 10 30 100" \
        --fo-values "1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1" \
        --formulation pressureJump --stab-factor 0.5 \
        --output-directory tauFo_pressureJump0p5 --clean
    ./undrained_mandel_sweep.py --tau-values "1 3 10 30 100" \
        --fo-values "1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1" \
        --formulation explicit --stab-factor 0.5 \
        --output-directory tauFo_explicit0p5 --clean

Each run writes a CSV summary and a four-panel map containing outer and
checkerboard contraction factors plus the checkerboard strength in the final
pressure field.  The preferred final-field diagnostic is the RMS of the
three-cell high-pass operator divided by pressure RMS; it rejects constant
and linear resolved trends while retaining grid-scale alternation.  The
preferred `interior_gridscale_to_pressure_rms` value additionally excludes
the outer two cell layers, so the prescribed top-boundary response is not
classified as bulk checkerboarding.  The CSV also records the corresponding
ratios to `pressure_fluctuation_rms`.  Those ratios expose contamination of
the nonuniform pressure response, while the pressure-RMS ratios measure the
actual perturbation relative to the full undrained pressure, including its
mean.  No pressure reference is used: all
sampled storage values are strictly positive and determine the pressure mean.

Matched tau--Fo directories can be compared directly without solver reruns:

    ./compare_results.py tauFo_unstabilized tauFo_stab0p5 \
        tauFo_explicit0p5 tauFo_pressureJump0p5 \
        --output-directory tauFoComparison

Use `--reuse-existing` on the sweep runner after changing diagnostics or plot
formatting.  Use `--rerun-unresolved` with a larger `--max-iterations` to
retain converged cases and repeat only points that previously reached the
iteration cap.

Stabilization-factor directories at one matched point can be summarized with
`plot_factor_sensitivity.py`; the CSV records both `formulation` and
`stab_factor` for provenance.
