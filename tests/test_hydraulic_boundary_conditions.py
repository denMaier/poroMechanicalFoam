from __future__ import annotations

import unittest

from tests.support import ROOT


class HydraulicBoundaryConditionTests(unittest.TestCase):
    def test_fixed_poro_flux_exposes_and_uses_minimum_conductivity_floor(self) -> None:
        header = (
            ROOT
            / "poroFluidModels/poroHydraulicFvPatchFields/fixedPoroFlux/fixedPoroFluxFvPatchScalarField.H"
        ).read_text(errors="ignore")
        source = (
            ROOT
            / "poroFluidModels/poroHydraulicFvPatchFields/fixedPoroFlux/fixedPoroFluxFvPatchScalarField.C"
        ).read_text(errors="ignore")

        self.assertIn("scalar minKeff_;", header)
        self.assertIn('lookupOrDefault<scalar>("minKeff", SMALL)', source)
        self.assertEqual(2, source.count("max(k_eff_, minKeff_)"))
        self.assertIn('writeKeyword("minKeff")', source)


if __name__ == "__main__":
    unittest.main()
