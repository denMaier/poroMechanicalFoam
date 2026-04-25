from __future__ import annotations

import unittest

from tests.support import ROOT, active_makefile_sources, strip_c_comments


class BuildManifestTests(unittest.TestCase):
    def test_active_sources_in_makefiles_exist(self) -> None:
        makefiles = [
            ROOT / "Make/files",
            ROOT / "materialModels/poroHydraulicModel/Make/files",
            *sorted((ROOT / "solvers").glob("*/Make/files")),
            *sorted((ROOT / "solvers/utilities").glob("*/Make/files")),
        ]

        missing: list[str] = []
        for makefile in makefiles:
            for source in active_makefile_sources(makefile):
                if not source.exists():
                    missing.append(f"{makefile.relative_to(ROOT)} -> {source.relative_to(ROOT)}")

        self.assertEqual([], missing)

    def test_makefiles_do_not_build_template_placeholders(self) -> None:
        offenders: list[str] = []
        for makefile in ROOT.rglob("Make/files"):
            active_text = strip_c_comments(makefile.read_text())
            if "[template]" in active_text or "templatePoroFluid" in active_text:
                offenders.append(str(makefile.relative_to(ROOT)))

        self.assertEqual([], offenders)

    def test_default_solver_build_matches_documented_targets(self) -> None:
        allwmake = (ROOT / "solvers/Allwmake").read_text()

        self.assertIn("cd poroMechanicalFoam", allwmake)
        self.assertIn("cd ../initPoroMechanicalFoam", allwmake)
        self.assertNotIn("poroMechanicalStrengthReduction && wmake", allwmake)

    def test_top_level_build_rejects_unsupported_openfoam_versions(self) -> None:
        allwmake = (ROOT / "Allwmake").read_text()

        self.assertIn('WM_PROJECT_VERSION" != "v2412"', allwmake)
        self.assertIn('WM_PROJECT_VERSION" != "v2512"', allwmake)
        self.assertIn("Unsupported OpenFOAM version detected", allwmake)


if __name__ == "__main__":
    unittest.main()
