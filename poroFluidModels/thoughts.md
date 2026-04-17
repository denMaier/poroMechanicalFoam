# thoughts

These are robustness improvements I noticed while documenting the `poroFluidModels` hierarchy. I did not change behavior here because the request was documentation-focused.

## 1. Guard assumptions around solver-performance lookup

Both variably saturated solvers read:

- `mesh().data().solverPerformanceDict().findEntry("p_rgh")`
- `mesh().data().solverPerformanceDict().findEntry("pHead")`

That assumes the dictionary entry always exists after the solve. A more robust version would check for the entry and fall back to a warning instead of dereferencing unconditionally.

## 2. Protect `newDeltaT()` against zero or degenerate coefficients

The adaptive time-step estimate divides by effective diffusivity-like terms such as:

- `Ss_ * magGamma`
- `(Ss_ + C) * magGamma`

If a material law or initialization produces very small or zero values, the estimate can become unstable or produce extreme time-step jumps. Clamping denominators and checking positivity explicitly would make this safer.

## 3. Make lazy initialization consistent

`poroFluid::poroHydraulic()` throws a fatal error if the pointer was not initialized earlier, while the variably saturated branch lazily creates its constitutive and linearization objects on first access.

Making the saturated branch follow the same lazy-init pattern would reduce ordering sensitivity and make the hierarchy more uniform.

## 4. Reduce dependence on string-matched residual names

`checkMassBalance()` only creates the residual field when iteration control already contains an entry named `"MassBalance"`.

That coupling-by-string is brittle. A more explicit switch in the dictionary, or a direct query through the iteration-control API, would make the behavior easier to reason about and less sensitive to naming drift.

## 5. Tighten error messages

At least one fatal error path in `poroFluid.C` reports the wrong function/class name (`varSatPoroFluidHead::poroHydraulic()`).

Cleaning those up would make failure reports more trustworthy when debugging initialization problems.

## 6. Normalize transient helper-field names

There are a few internal names that look accidental or debug-oriented, for example the `varSatPoroFluidHead` helper field is created with the IO name `"makeKEfff"`.

Using stable semantic field names for internal helper objects would make debugging and restart behavior less surprising.
