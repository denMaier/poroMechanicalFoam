# Fully undrained strip-footing stability test

This test is the mechanically forced counterpart of `undrainedMandel`.  It
uses the same material, mesh spacing, coupling diagnostics, and definitions

    tau = M/K_vol,
    Fo  = Kp dt / ((Ss + 1/K_vol) h^2),
    Ss  = 1/(tau K_vol).

The domain represents one symmetric half of a strip footing.  The boundary
`x=0` is a symmetry plane.  A total compressive traction is prescribed over
`0 <= x <= 0.25` on `topLoaded`; the remainder of the top is traction-free.
The bottom is fixed and the right side is traction-free.  Every hydraulic
boundary has zero normal flux, and `Ss` is always positive, so no pressure
reference is used.

The sharp change from loaded to unloaded top surface intentionally reproduces
the standard footing test used to expose pressure instability in nearly
impermeable poroelasticity.  It is a mechanical traction discontinuity rather
than a prescribed pressure or flux boundary layer.  A smoothed-load variant
should subsequently be used to check sensitivity to this load edge.

Run an unstabilized pilot in the low-regularization corner with

    ./undrained_footing_sweep.py --tau-values "100" --fo-values "1e-6" \
        --formulation unstabilized --max-iterations 1200 \
        --output-directory pilot_unstabilized --clean

The default footing traction is `-1e5 Pa`; it can be changed using
`--footing-traction`.  Full matched sweeps use the same commands as the Mandel
test, replacing the runner and output-directory names.

The CSV reports total, pressure-equation, and checkerboard contraction rates.
Final pressure checkerboarding is measured both relative to total pressure RMS
and relative to the smooth pressure fluctuation RMS.  The interior high-pass
metric excludes two cell layers at the exterior boundaries.
