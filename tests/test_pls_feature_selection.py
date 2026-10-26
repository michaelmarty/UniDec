import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd


class TestPLSFeatureSelection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = (Path(__file__).resolve().parents[3] / "Scripts/Annika"
                / "DDA Paper Scripts/PLS_Feature_Selection.py")
        spec = importlib.util.spec_from_file_location("pls_selection", path)
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def test_grouped_folds_keep_repeat_lipids_together(self):
        groups = np.repeat(np.arange(100), 2)
        y = np.repeat(np.arange(100) % 2, 2)
        seen = np.zeros(len(y))
        for train, test in self.module.grouped_folds(y, groups):
            self.assertFalse(set(groups[train]) & set(groups[test]))
            seen[test] += 1
        np.testing.assert_array_equal(seen, np.ones(len(y)))

    def test_all_triplets_use_identical_held_out_rows(self):
        generator = np.random.default_rng(42)
        groups = np.arange(100)
        y = groups % 2
        x = np.column_stack([y + generator.normal(0, 0.1, 100),
                             generator.normal(size=(100, 3))])
        best, ranking, predictions = self.module.rank_triplets(
            x, y, ["signal", "noise1", "noise2", "noise3"],
            self.module.grouped_folds(y, groups))
        self.assertEqual(len(ranking), 4)
        self.assertIn("signal", best)
        for prediction in predictions.values():
            self.assertTrue(np.isfinite(prediction).all())
        self.assertGreater(ranking.iloc[0]["CV AUC"], 0.99)

    def test_candidate_pool_excludes_status_and_uses_common_rows(self):
        frame = pd.DataFrame({
            "Ontology": ["PC"] * 6, "Metabolite name": [f"PC {i}" for i in range(6)],
            "Overall Match Status": [True, False, True, False, "Provisional", True],
            "Square Root Class Strength": [1, 2, 3, 4, 5, np.nan],
            "Reverse dot product": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            "Fill %": [1, 0.33, 1, 0.67, 1, 1],
            "Matched Required Fragments": ["Yes"] * 6,
            "Matched peaks count": [1, 2, 3, 4, 5, 6],
            "Matched peaks percentage": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        })
        rows, names, _ = self.module.prepare_candidates(frame)
        self.assertEqual(len(rows), 4)
        self.assertNotIn("Overall Match Status", names)
        self.assertNotIn("Matched Required Fragments", names)
        self.assertNotIn("Matched peaks count", names)
        self.assertNotIn("Matched peaks percentage", names)


if __name__ == "__main__":
    unittest.main()
