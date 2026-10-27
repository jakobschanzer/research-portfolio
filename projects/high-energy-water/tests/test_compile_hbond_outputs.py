from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "python"
    / "compile_hbond_outputs.py"
)
SPEC = importlib.util.spec_from_file_location("compile_hbond_outputs", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CompileHydrogenBondOutputsTest(unittest.TestCase):
    def write_fixture(self, directory: Path) -> None:
        for temperature in MODULE.EXPECTED_TEMPERATURES:
            values = [
                f"{step} {temperature / 10 + step:.2f} {step + 0.25:.2f} {step + 0.75:.2f}"
                for step in range(1, 5)
            ]
            (directory / f"dist.{temperature}").write_text(
                "step HB_WW HB_WH_top HB_WH_bottom\n" + "\n".join(values) + "\n",
                encoding="utf-8",
            )

    def test_seeded_compilation_is_reproducible(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs = root / "inputs"
            first = root / "first"
            second = root / "second"
            inputs.mkdir()
            self.write_fixture(inputs)

            outputs_a = MODULE.compile_outputs(inputs, first, seed=17, iterations=100)
            outputs_b = MODULE.compile_outputs(inputs, second, seed=17, iterations=100)

            for path_a, path_b in zip(outputs_a, outputs_b, strict=True):
                self.assertEqual(path_a.read_bytes(), path_b.read_bytes())

            with outputs_a[1].open(encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 10)
            self.assertEqual(rows[0]["Temp_C"], "5.0")
            self.assertAlmostEqual(float(rows[0]["HB_WW"]), 3.0)

    def test_missing_temperature_is_reported(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaisesRegex(FileNotFoundError, "missing temperature files"):
                MODULE.compile_outputs(root, root / "output", iterations=10)


if __name__ == "__main__":
    unittest.main()

