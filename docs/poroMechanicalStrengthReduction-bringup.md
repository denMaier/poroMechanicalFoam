# Bringing Up `poroMechanicalStrengthReduction`

This note describes the work required to get `poroMechanicalStrengthReduction`
to the point where it can be compiled, run on a prepared case, and then
validated as a usable strength reduction workflow.

## Current status

`solvers/poroMechanicalStrengthReduction/` already contains a real solver
implementation, not just a stub:

- `poroMechanicalStrengthReduction.C` reads `constant/mechanicalProperties`
- it expects `startFrictionAngle` and `startCohesion`
- it updates `frictionAngle`, `cohesion`, and `dilationAngle`
- it creates a `physicsModel` and advances the case step-by-step

However, it is not yet production-ready:

- the source file is marked `UNFINISHED WORK!`
- the solver is not built by `solvers/Allwmake`
- no bundled tutorial or example case for this solver was found
- no recent build log for this solver was found in the repository

## Goal levels

There are two distinct bring-up targets:

1. `Compile + first run`
   Get the solver to build and execute on one curated case.
2. `Usable workflow`
   Make it reliable, documented, and validated enough for repeated use.

The first target looks feasible with moderate work. The second target will
require additional implementation and validation.

## Required work for `compile + first run`

### 1. Build integration

The solver must be added to the normal solver build path.

Required actions:

- update `solvers/Allwmake` so it also builds `poroMechanicalStrengthReduction`
- confirm `Make/files` symlink handling is correct for both OpenFOAM and
  foam-extend modes
- verify the executable lands in `$(FOAM_USER_APPBIN)`

Expected result:

- `./Allwmake` builds `poroMechanicalStrengthReduction` together with the other
  solvers

### 2. Compile verification against supported OpenFOAM versions

The repository README states support for OpenFOAM `v2412` and `v2512`. The
solver must be compiled against at least one supported environment first, and
preferably both.

Required actions:

- source a supported OpenFOAM environment
- run the repository `./Allwmake`
- capture any compile errors from `solvers/poroMechanicalStrengthReduction`
- fix API drift or missing includes if they appear

Likely focus areas:

- file copy helper usage (`cp(...)`)
- `physicsModel` lifecycle assumptions
- dictionary write/update behavior for `mechanicalProperties`
- compatibility of linked libraries in
  `solvers/poroMechanicalStrengthReduction/Make/options`

Expected result:

- the solver compiles cleanly in at least one supported OpenFOAM version

### 3. Minimal runnable case

No bundled tutorial case for this solver was found, so a case must be prepared
 explicitly for it.

Required actions:

- choose a working `poroMechanicalFoam` case as the starting point
- switch `application` in `system/controlDict` to
  `poroMechanicalStrengthReduction`
- ensure `constant/mechanicalProperties` contains a compatible mechanical law
- add the required entries:
  - `startFrictionAngle`
  - `startCohesion`
- confirm the solver is allowed to update:
  - `frictionAngle`
  - `cohesion`
  - `dilationAngle`

Important assumption in current code:

- the solver reads the first entry of the `mechanical` list and modifies that
  dictionary directly

This means the selected case must match that expectation.

Expected result:

- one case can be launched without immediate dictionary or field lookup errors

### 4. Field continuity between time steps

The solver currently copies `TotalHead` and `Sw` from the previous time
directory into the next one after each solve step.

Required actions:

- confirm the target case actually writes `TotalHead`
- confirm the target case actually writes `Sw`
- verify both files exist at every step where the copy is attempted
- decide what should happen for fully saturated or differently configured cases

Potential implementation work:

- guard the copy operations if either field is absent
- generalize the restart/field transfer logic if more fields are required

Expected result:

- the solver can advance multiple time steps without failing on missing field
  files

## Work required for a usable strength reduction workflow

The current solver appears to reduce parameters using time progression as a
factor-of-safety sequence, but a reliable workflow needs more than that.

### 5. Failure criterion

This is the main functional gap.

The current implementation reduces strength parameters over steps, but there is
no obvious automated failure detection or stopping criterion in the solver.

Required decisions:

- define what constitutes failure
- define when the run should stop
- define what FoS value should be reported

Possible criteria:

- non-convergence of the mechanical solve
- divergence of displacement or residual metrics
- second-order work or other instability indicator
- loss of equilibrium over a prescribed tolerance

Expected result:

- the solver reports a meaningful failure FoS instead of only marching through
  time-like steps

### 6. Output and reporting

A strength reduction solver should make the result easy to inspect.

Required actions:

- define what scalar result is the main output
- write the current reduction state clearly at each step
- write a final summary with:
  - last converged FoS
  - failed FoS
  - selected failure criterion

Expected result:

- a run produces a clear, defensible FoS outcome

### 7. Case documentation

Without a tutorial or recipe, the solver is difficult to use and difficult to
maintain.

Required actions:

- add one documented example case
- document required entries in `mechanicalProperties`
- document expected `controlDict` stepping
- document limitations and assumptions

Expected result:

- another user can reproduce the workflow without reading solver source first

### 8. Validation

A compile-successful solver is not yet trustworthy.

Required actions:

- compare results against a known strength reduction benchmark
- verify sensitivity to:
  - time-step / FoS increment size
  - initial stress state
  - hydraulic field initialization
  - constitutive law choice
- check whether copying `TotalHead` and `Sw` between steps is physically and
  numerically consistent for the intended workflow

Expected result:

- confidence that the reported FoS is not just a numerical artifact

## Suggested implementation order

1. Add the solver to `solvers/Allwmake`
2. Compile it in a sourced OpenFOAM environment
3. Fix any compile errors
4. Prepare one minimal case with the required mechanical entries
5. Run until the first successful multi-step execution
6. Harden field-copy and dictionary assumptions
7. Implement or formalize failure detection
8. Add a documented validation case

## Concrete unknowns to resolve during bring-up

These are the items most likely to require actual coding rather than just setup:

- whether `cp(...)` is valid and portable in the targeted OpenFOAM version(s)
- whether rewriting `mechanicalProperties` at runtime is robust across all
  cases
- whether `physicsModel::New(runTime)` each loop iteration is intended and safe
- whether all necessary restart fields are covered by only copying `TotalHead`
  and `Sw`
- whether the solver should always write every step or only output steps
- how failure is detected and how the final FoS is selected

## Practical conclusion

Bringing `poroMechanicalStrengthReduction` to the point of compiling and
running on a prepared case looks feasible without a major rewrite.

Bringing it to the point of being a reliable and documented strength reduction
tool will require additional implementation work, mainly around:

- build integration
- case setup definition
- restart/field handling
- failure detection
- validation
