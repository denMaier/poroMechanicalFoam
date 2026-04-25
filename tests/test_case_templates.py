from __future__ import annotations

import unittest

from tests.support import ROOT, dictionary_value, header_typenames


CASE_ROOT = ROOT / "caseTemplates"


class CaseTemplateTests(unittest.TestCase):
    def test_every_case_template_has_minimum_openfoam_layout(self) -> None:
        missing: list[str] = []

        for case in sorted(path for path in CASE_ROOT.iterdir() if path.is_dir()):
            required = [
                case / "0",
                case / "constant",
                case / "system/controlDict",
                case / "system/fvSchemes",
                case / "system/fvSolution",
                case / "constant/g",
                case / "constant/poroFluidProperties",
                case / "constant/poroHydraulicProperties",
            ]

            if "poroMechanical" in case.name:
                required.extend(
                    [
                        case / "constant/mechanicalProperties",
                        case / "constant/poroCouplingProperties",
                    ]
                )
            else:
                required.append(case / "constant/polyMesh/blockMeshDict")

            for path in required:
                if not path.exists():
                    missing.append(str(path.relative_to(ROOT)))

        self.assertEqual([], missing)

    def test_case_templates_use_the_main_solver(self) -> None:
        bad_applications: list[str] = []

        for control_dict in sorted(CASE_ROOT.glob("*/system/controlDict")):
            application = dictionary_value(control_dict, "application")
            if application != "poroMechanicalFoam":
                bad_applications.append(
                    f"{control_dict.relative_to(ROOT)} declares {application!r}"
                )

        self.assertEqual([], bad_applications)

    def test_template_runtime_model_names_are_registered_in_repo(self) -> None:
        known_types = header_typenames()
        expected_keys = {
            "poroFluidModel",
            "poroSolidInterface",
            "storageLaw",
            "saturationLaw",
            "effectiveStressModel",
        }

        missing: list[str] = []
        for dictionary in sorted(CASE_ROOT.glob("*/constant/*Properties")):
            for key in expected_keys:
                model_name = dictionary_value(dictionary, key)
                if model_name and model_name not in known_types:
                    missing.append(f"{dictionary.relative_to(ROOT)}: {key} {model_name}")

        self.assertEqual([], missing)

    def test_flow_templates_have_the_field_selected_by_their_model(self) -> None:
        expectations = {
            "poroFluid-saturated": "p_rgh",
            "poroFluid-varSat-p_rgh": "p_rgh",
            "poroFluid-varSat-pHead": "pHead",
        }

        missing: list[str] = []
        for case_name, field_name in expectations.items():
            field = CASE_ROOT / case_name / "0" / field_name
            if not field.exists():
                missing.append(str(field.relative_to(ROOT)))

        self.assertEqual([], missing)


if __name__ == "__main__":
    unittest.main()
