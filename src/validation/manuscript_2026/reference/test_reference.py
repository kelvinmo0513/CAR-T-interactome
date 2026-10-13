"""Small synthetic fixtures; no real cell identifiers or clinical records."""
import unittest
from pathlib import Path

from oversample import TARGET, expansion, require_external_output


class ReferenceExpansionTests(unittest.TestCase):
    def fixture(self):
        return [dict(source_cell="synthetic_a", source_group="synthetic_sample", phenotype="B cells"),
                dict(source_cell="synthetic_b", source_group="synthetic_sample", phenotype=TARGET),
                dict(source_cell="synthetic_c", source_group="synthetic_sample", phenotype="T cells"),
                dict(source_cell="synthetic_d", source_group="synthetic_sample", phenotype=TARGET)]

    def test_target_expanded_in_place(self):
        self.assertEqual(expansion(self.fixture(), 4), [0] + [1] * 10 + [2] + [3] * 10)

    def test_duplicate_source_id_rejected(self):
        fixture = self.fixture()
        fixture[3]["source_cell"] = fixture[1]["source_cell"]
        with self.assertRaisesRegex(ValueError, "unique"):
            expansion(fixture, 4)

    def test_single_target_rejected(self):
        fixture = self.fixture()
        fixture[3]["phenotype"] = "Tregs"
        with self.assertRaisesRegex(ValueError, "Exactly two"):
            expansion(fixture, 4)

    def test_changed_source_support_rejected(self):
        fixture = self.fixture()
        fixture[3]["source_group"] = "different_synthetic_sample"
        with self.assertRaisesRegex(ValueError, "support differs"):
            expansion(fixture, 4)

    def test_repository_output_rejected(self):
        with self.assertRaisesRegex(ValueError, "outside the repository"):
            require_external_output(Path(__file__).parent / "private_output")


if __name__ == "__main__":
    unittest.main()
