from __future__ import annotations

import unittest

from tests.support import ROOT


class PoroTractionConsistencyTests(unittest.TestCase):
    def test_mechanical_law_exposes_its_pressure_field(self) -> None:
        header = (
            ROOT
            / "materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/"
            "poroMechanicalLaw2/poroMechanicalLaw2.H"
        ).read_text(errors="ignore")

        self.assertIn("const word& pressureFieldName() const", header)
        self.assertIn("return pName_;", header)
        self.assertIn("const dimensionedScalar& biotCoeffValue() const", header)

    def test_effective_traction_uses_biot_weighted_law_pressure(self) -> None:
        source = (
            ROOT
            / "poroSolidInterfaces/poroSolidFvPatchFields/poroTraction/"
            "poroTractionFvPatchVectorField.C"
        ).read_text(errors="ignore")

        self.assertIn('pressureFieldName_(dict.lookupOrDefault<word>("pressureFieldName", "auto"))', source)
        self.assertIn("const word pName(mechanicalPressureFieldName(solMod));", source)
        self.assertIn("const scalarField bPatch(patchBiotCoeff(solMod));", source)
        self.assertIn("totalPressure = pressure() + bPatch*pPatch;", source)

    def test_explicit_pressure_selection_is_checked_against_the_law(self) -> None:
        source = (
            ROOT
            / "poroSolidInterfaces/poroSolidFvPatchFields/poroTraction/"
            "poroTractionFvPatchVectorField.C"
        ).read_text(errors="ignore")

        self.assertIn("pressureFieldName_ != lawPressureName", source)
        self.assertIn("Pressure-field mismatch on poroTraction patch", source)
        self.assertIn("Entry 'buoyancyIncluded'", source)
        self.assertIn("is deprecated", source)


if __name__ == "__main__":
    unittest.main()
