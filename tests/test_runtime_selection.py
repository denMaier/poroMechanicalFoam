from __future__ import annotations

import re
import unittest

from tests.support import ROOT, active_makefile_sources, header_typenames


REGISTRATION_RE = re.compile(
    r"addToRunTimeSelectionTable\s*\(\s*([A-Za-z_]\w*)\s*,\s*([A-Za-z_]\w*)",
    flags=re.DOTALL,
)


class RuntimeSelectionTests(unittest.TestCase):
    def test_compiled_runtime_selected_classes_define_type_names(self) -> None:
        makefiles = [
            ROOT / "Make/files",
            ROOT / "materialModels/poroHydraulicModel/Make/files",
        ]

        missing_definitions: list[str] = []
        for makefile in makefiles:
            for source in active_makefile_sources(makefile):
                text = source.read_text(errors="ignore")
                for _, class_name in REGISTRATION_RE.findall(text):
                    if f"defineTypeNameAndDebug({class_name}" not in text:
                        missing_definitions.append(
                            f"{source.relative_to(ROOT)} registers {class_name}"
                        )

        self.assertEqual([], missing_definitions)

    def test_runtime_selected_classes_have_matching_header_typename(self) -> None:
        missing_headers: list[str] = []
        makefiles = [
            ROOT / "Make/files",
            ROOT / "materialModels/poroHydraulicModel/Make/files",
        ]

        for makefile in makefiles:
            for source in active_makefile_sources(makefile):
                text = source.read_text(errors="ignore")
                for _, class_name in REGISTRATION_RE.findall(text):
                    headers = list(source.parent.glob("*.H"))
                    matching_headers = [
                        path
                        for path in headers
                        if re.search(
                            rf"\bclass\s+{re.escape(class_name)}\b",
                            path.read_text(errors="ignore"),
                        )
                    ]

                    if not matching_headers:
                        missing_headers.append(
                            f"{source.relative_to(ROOT)} registers {class_name} without a neighboring class header"
                        )
                        continue

                    if not any(
                        "TypeName(" in path.read_text(errors="ignore")
                        for path in matching_headers
                    ):
                        locations = ", ".join(
                            str(path.relative_to(ROOT)) for path in matching_headers
                        )
                        missing_headers.append(
                            f"{source.relative_to(ROOT)} registers {class_name}; no TypeName in {locations}"
                        )

        self.assertEqual([], missing_headers)

    def test_richards_linearization_default_build_excludes_newton_consistently(self) -> None:
        makefile = ROOT / "Make/files"
        active_sources = {path.relative_to(ROOT).as_posix() for path in active_makefile_sources(makefile)}

        self.assertIn(
            "poroFluidModels/richardsLinearization/Standard/Standard.C",
            active_sources,
        )
        self.assertNotIn(
            "poroFluidModels/richardsLinearization/Newton/Newton.C",
            active_sources,
        )

    def test_registered_type_names_are_unique_within_model_family_directories(self) -> None:
        family_roots = [
            ROOT / "poroFluidModels",
            ROOT / "poroSolidInterfaces",
            ROOT / "materialModels/poroHydraulicModel/saturationLaws",
            ROOT / "materialModels/poroHydraulicModel/storageLaws",
            ROOT / "materialModels/poroHydraulicModel/conductivityModels",
            ROOT / "materialModels/mechanicalModel/mechanicalLaws/linearGeometryLaws/varSatPoroMechanicalLaw/effectiveStressModels",
            ROOT / "iterationControls",
        ]

        duplicates: list[str] = []
        for family_root in family_roots:
            typenames = header_typenames(list(family_root.rglob("*.H")))
            for type_name, paths in sorted(typenames.items()):
                if len(paths) > 1:
                    location_list = ", ".join(str(path.relative_to(ROOT)) for path in sorted(paths))
                    duplicates.append(f"{family_root.relative_to(ROOT)}: {type_name} in {location_list}")

        self.assertEqual([], duplicates)


if __name__ == "__main__":
    unittest.main()
