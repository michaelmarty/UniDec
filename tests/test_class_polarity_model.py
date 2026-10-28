import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd


class TestClassPolarityModel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = (Path(__file__).resolve().parents[3] / "Scripts/Annika"
                / "DDA Paper Scripts/Class_Polarity_Model.py")
        spec = importlib.util.spec_from_file_location("class_polarity_model", path)
        cls.module = importlib.util.module_from_spec(spec)
        with patch.object(sys, "path", [str(path.parent)] + sys.path):
            spec.loader.exec_module(cls.module)

    def fixture(self):
        return pd.DataFrame([
            {"Ontology": "PC", "Dataset": "Stellar", "Polarity": mode,
             "Metabolite name": f"PC {identity}:0", "Average Rt(min)": identity / 10 + offset,
             "Overall Match Status": bool(identity % 2),
             "Class Strength": (9 if identity % 2 else 1) if mode == "Positive"
             else (1 if identity % 2 else 9)}
            for identity in range(80)
            for mode, offset in [("Positive", 0), ("Negative", 0.05)]
        ])

    def test_max_normalization_and_zero_group(self):
        rows = pd.DataFrame({"Ontology": ["PC", " pc ", "PC", "PE", "PE"],
                             "Polarity": ["Positive", "positive", "Negative", "Negative", "Negative"],
                             "Class Strength": [2, 4, 10, 0, np.nan]})
        result = self.module.normalize_strength(rows)
        np.testing.assert_allclose(result[self.module.NORMALIZED], [0.5, 1, 1, 0, np.nan])

    def test_training_maximum_excludes_test_scores(self):
        frame = self.module.prepare_rows(self.fixture())
        frame.loc[0, "Class Strength"] = 1000
        score, _ = self.module.scaled_scores(frame, np.arange(2, len(frame)))
        self.assertAlmostEqual(score[0], 1000 / 9)

    def test_matching_rejects_ties_and_cross_dataset_pairs(self):
        rows = self.fixture().iloc[:4].copy()
        duplicate = rows.iloc[[1]].copy()
        rows.loc[3, "Dataset"] = "OTOT"
        frame = self.module.prepare_rows(pd.concat([rows, duplicate], ignore_index=True))
        self.assertTrue(frame["Partner Row"].eq(-1).all())

    def test_group_specific_direction_and_held_out_identity(self):
        frame = self.module.prepare_rows(self.fixture())
        summary, predictions, cutoffs = self.module.evaluate(frame)
        aucs = summary.set_index("Model")["CV AUC"]
        self.assertGreater(aucs["Class/polarity score"], 0.99)
        self.assertGreater(aucs["Combined polarities"], 0.99)
        self.assertLess(aucs["Global score"], 0.6)
        self.assertEqual(predictions.groupby("Lipid Group")["Fold"].nunique().max(), 1)
        self.assertEqual(set(cutoffs["Direction"]), {"Above", "Below"})
        self.assertTrue(frame["Partner Row"].ge(0).all())


if __name__ == "__main__":
    unittest.main()
