import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd


class TestClassStrengthAnalysis(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = (Path(__file__).resolve().parents[3] / "Scripts/Annika"
                / "DDA Paper Scripts/Class_Strength_Analysis.py")
        spec = importlib.util.spec_from_file_location("strength_analysis", path)
        cls.module = importlib.util.module_from_spec(spec)
        with patch.object(sys, "path", [str(path.parent)] + sys.path):
            spec.loader.exec_module(cls.module)

    def fixture(self):
        return pd.DataFrame([
            {"Ontology": "PC", "Polarity": "Positive", "Metabolite name": f"PC {identity}",
             "Overall Match Status": identity % 2 == 1,
             "Square Root Class Strength": 10 if identity % 2 else 1}
            for identity in range(80) for _ in range(2)
        ])

    def test_strength_only_model_needs_no_other_metric(self):
        frame, model, auc = self.module.fit_class_strength(self.fixture())
        self.assertEqual(len(frame), 160)
        self.assertEqual(model.n_features_in_, 1)
        self.assertGreater(auc, 0.99)
        probabilities = frame["Class Strength OOF Probability"]
        self.assertTrue(probabilities.between(0, 1).all())
        self.assertEqual(frame.groupby("Lipid Group")["Fold"].nunique().max(), 1)

    def test_provisional_and_nonfinite_strength_are_excluded(self):
        rows = self.fixture()
        rows["Overall Match Status"] = rows["Overall Match Status"].astype(object)
        rows["Square Root Class Strength"] = rows["Square Root Class Strength"].astype(float)
        rows.loc[0, "Overall Match Status"] = "Provisional"
        rows.loc[1, "Square Root Class Strength"] = np.inf
        rows.loc[2, "Square Root Class Strength"] = np.nan
        frame, _, _ = self.module.fit_class_strength(rows)
        self.assertEqual(len(frame), 157)
        self.assertTrue(np.isfinite(frame["Square Root Class Strength"]).all())

    def test_benchmark_compares_same_rows_without_outcome_predictors(self):
        rows = self.fixture()
        rng = np.random.default_rng(42)
        rows["Reverse dot product"] = rng.random(len(rows))
        rows["Fill %"] = rng.random(len(rows))
        rows.loc[0, "Fill %"] = np.nan
        summary, predictions = self.module.benchmark_metrics(rows)
        self.assertEqual(set(summary["Rows"]), {159})
        self.assertEqual(len(predictions), 159)
        self.assertEqual(summary.iloc[0]["Metric"], self.module.STRENGTH)
        self.assertNotIn("Overall Match Status", set(summary["Metric"]))


if __name__ == "__main__":
    unittest.main()
