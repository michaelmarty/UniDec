import importlib.util
from pathlib import Path
import unittest

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


class TestSIClassScoreComparison(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        script = (Path(__file__).resolve().parents[3] / "Scripts" / "Annika"
                  / "DDA Paper Scripts" / "SI_Class_Score_Comparison.py")
        spec = importlib.util.spec_from_file_location("si_class_scores", script)
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def fixture(self):
        rows = []
        for identity in range(40):
            accepted = identity >= 20
            for replicate in range(2):
                rows.append({
                    "Ontology": "SM", "Polarity": "Negative",
                    "Metabolite name": f"SM {identity}:0|SM {identity}:0_{replicate}:0",
                    "Overall Match Status": accepted,
                    "Square Root Class Strength": 5,
                    "Reverse dot product": 0.1 if accepted else 0.9,
                    "Fill %": 80 if accepted else 20,
                })
        return pd.DataFrame(rows)

    def test_provisional_rows_are_not_negative_labels(self):
        rows = self.fixture()
        rows["Overall Match Status"] = rows["Overall Match Status"].astype(object)
        rows.loc[0, "Overall Match Status"] = "Provisional"
        rows.loc[1, "Overall Match Status"] = None
        rows["Validation Status"] = "Accepted"
        rows.loc[2, "Validation Status"] = "Provisional"
        frame = self.module.prepare_scores(rows, ["SM"])
        self.assertEqual(len(frame), len(rows) - 3)
        self.assertEqual(int((~frame["Accepted"]).sum()), 37)

    def test_repeated_identities_do_not_cross_train_test_boundary(self):
        frame = self.module.prepare_scores(self.fixture(), ["SM"])
        y = frame["Accepted"].to_numpy(dtype=int)
        groups = frame["Lipid Group"].to_numpy()
        seen = np.zeros(len(frame), dtype=int)
        for train, test in self.module.grouped_folds(y, groups):
            self.assertFalse(set(groups[train]) & set(groups[test]))
            self.assertEqual(len(np.unique(y[train])), 2)
            self.assertEqual(len(np.unique(y[test])), 2)
            seen[test] += 1
        np.testing.assert_array_equal(seen, np.ones(len(frame)))

    def test_model_learns_reversed_direction_on_training_data(self):
        frame = self.module.prepare_scores(self.fixture(), ["SM"])
        summary, predictions = self.module.evaluate_scores(frame)
        table = summary.set_index("Score")
        self.assertAlmostEqual(table.loc["Class strength", "CV AUC"], 0.5)
        self.assertAlmostEqual(table.loc["Reverse dot product", "CV AUC"], 1)
        self.assertAlmostEqual(table.loc["Fill %", "CV AUC"], 1)
        self.assertEqual(predictions.groupby("Lipid Group")["Fold"].nunique().max(), 1)
        self.assertTrue((table["CI low"] <= table["CI high"]).all())
        self.assertGreater(table.loc["Reverse dot product", "Delta CI low"], 0)

    def test_one_outcome_is_not_given_a_prediction_auc(self):
        source = self.fixture()
        source["Overall Match Status"] = True
        frame = self.module.prepare_scores(source, ["SM"])
        summary, predictions = self.module.evaluate_scores(frame)
        self.assertTrue(summary["CV AUC"].isna().all())
        self.assertTrue(summary["Folds"].eq(0).all())
        self.assertTrue(predictions.empty)

    def test_missing_score_uses_identical_rows_for_all_models(self):
        source = self.fixture()
        source.loc[0, "Fill %"] = np.nan
        frame = self.module.prepare_scores(source, ["SM"])
        self.assertEqual(len(frame), len(source) - 1)
        self.assertTrue(np.isfinite(frame[list(self.module.SCORES)]).all().all())


if __name__ == "__main__":
    unittest.main()
