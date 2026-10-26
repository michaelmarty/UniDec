import importlib.util
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


class TestSIOverlapSN(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        script = (Path(__file__).resolve().parents[3] / "Scripts" / "Annika"
                  / "DDA Analysis" / "SI_Overlap_SN.py")
        spec = importlib.util.spec_from_file_location("si_overlap_sn", script)
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def tearDown(self):
        self.module.plt.close("all")

    def fixture(self):
        rows = [(dataset, "PC 34:1|PC 16:0_18:1", True, (index + 1) * 10)
                for index, dataset in enumerate(self.module.DATASETS)]
        rows += [
            (" OTOT\u00a0", "PC 34:1|other detail", "TRUE", 30),
            ("OTOT", "unique", True, 5),
            ("OTIT", "unique", False, 500),
            ("Apex", "unique", True, 100),
            ("OTOT", "shared", True, 8),
            ("OTIT", "shared", True, "bad"),
            ("ITIT", "zero", True, 0),
            ("Stellar", "infinite", True, np.inf),
            ("Stellarboth", "negative", True, -2),
            ("OTOT", None, True, 50),
        ]
        return pd.DataFrame(rows, columns=["Dataset", "Metabolite name",
                                           "Overall Match Status", "S/N average"])

    def test_lipid_identity_acceptance_and_method_specific_medians(self):
        original = self.fixture()
        lipids = self.module.prepare_lipid_sn(original)
        shared = lipids.loc[lipids["Lipid"].eq("PC 34:1")].set_index("Dataset")
        self.assertEqual(shared.loc["OTOT", "S/N"], 20)
        self.assertEqual(shared.loc["OTOT", "Accepted Rows"], 2)
        self.assertEqual(shared.loc["ITIT", "S/N"], 30)
        unique = lipids.loc[lipids["Lipid"].eq("unique")]
        self.assertEqual(unique["Dataset"].tolist(), ["OTOT"])
        self.assertFalse(lipids.duplicated(["Dataset", "Lipid"]).any())
        self.assertNotIn("Lipid", original)

    def test_invalid_signal_keeps_accepted_counts(self):
        lipids = self.module.prepare_lipid_sn(self.fixture())
        self.assertEqual(int(lipids["S/N"].isna().sum()), 4)
        summary = self.module.summarize_lipid_sn(lipids)
        self.assertEqual(len(summary), 5)
        self.assertEqual(summary.loc["OTIT", "Accepted_Lipids"], 2)
        self.assertEqual(summary.loc["OTIT", "Plotted_Lipids"], 1)
        self.assertEqual(summary.loc["Stellar", "Accepted_Lipids"], 2)

    def test_one_log_axis_has_one_dot_per_usable_lipid(self):
        lipids = self.module.prepare_lipid_sn(self.fixture())
        fig = self.module.plot_dataset_sn(lipids)
        self.assertEqual(len(fig.axes), 1)
        ax = fig.axes[0]
        self.assertEqual(ax.get_yscale(), "log")
        self.assertEqual([tick.get_text() for tick in ax.get_xticklabels()],
                         [self.module.DATASET_LABELS[name] for name in self.module.DATASETS])
        for dataset, collection in zip(self.module.DATASETS, ax.collections):
            plotted = collection.get_offsets()[:, 1]
            expected = lipids.loc[lipids["Dataset"].eq(dataset), "S/N"].dropna()
            np.testing.assert_allclose(sorted(plotted), sorted(expected))

    def test_empty_and_unusable_populations(self):
        frame = self.fixture()
        with self.assertRaisesRegex(ValueError, "Missing required columns"):
            self.module.prepare_lipid_sn(frame.drop(columns="S/N average"))
        with self.assertRaisesRegex(ValueError, "No accepted lipids"):
            self.module.prepare_lipid_sn(frame.iloc[:0])
        frame["S/N average"] = 0
        with self.assertRaisesRegex(ValueError, "No positive finite values"):
            self.module.plot_dataset_sn(self.module.prepare_lipid_sn(frame))

    def test_repeated_measures_anova_and_paired_comparisons(self):
        rng = np.random.default_rng(42)
        base = rng.normal(2, 0.5, 60)
        rows = [(dataset, str(lipid), 10 ** (base[lipid] + method * 0.3 + rng.normal(0, 0.15)))
                for method, dataset in enumerate(self.module.DATASETS)
                for lipid in range(60)]
        lipids = pd.DataFrame(rows, columns=["Dataset", "Lipid", "S/N"])
        anova, pairs = self.module.compare_dataset_sn(lipids)
        self.assertEqual(anova["df_num"], 4)
        self.assertEqual(anova["df_den"], 236)
        self.assertEqual(anova["Matched lipids"], 60)
        self.assertGreaterEqual(anova["GG epsilon"], 0.25)
        self.assertLessEqual(anova["GG epsilon"], 1)
        self.assertGreaterEqual(anova["p"], anova["p_uncorrected"])
        values = np.log10(lipids.pivot(index="Lipid", columns="Dataset", values="S/N").to_numpy())
        residual = values - values.mean(axis=0) - values.mean(axis=1)[:, None] + values.mean()
        expected_f = (60 * np.sum((values.mean(axis=0) - values.mean()) ** 2) / 4
                      / (np.sum(residual ** 2) / 236))
        np.testing.assert_allclose(anova["F"], expected_f)
        self.assertLess(anova["p"], 0.001)
        self.assertEqual(len(pairs), 10)
        self.assertTrue(pairs["Holm p"].ge(pairs["p"]).all())
        self.assertTrue(pairs["Significant"].all())
        np.testing.assert_allclose(pairs.iloc[0]["Geometric mean ratio"],
                                   10 ** (np.log10(lipids.loc[lipids.Dataset.eq("OTOT"), "S/N"]).mean()
                                          - np.log10(lipids.loc[lipids.Dataset.eq("OTIT"), "S/N"]).mean()))
        fig = self.module.plot_dataset_sn(lipids, (anova, pairs))
        self.assertTrue(all(ax.get_title() == "" for ax in fig.axes))
        self.assertEqual(len(fig.axes[1].lines), 0)
        self.assertFalse(any(text.get_text().startswith("n=") for text in fig.axes[0].texts))
        with self.assertRaisesRegex(ValueError, "at least two lipids"):
            self.module.compare_dataset_sn(lipids.loc[lipids.Dataset.ne("OTOT")])
        with self.assertRaisesRegex(ValueError, "one value per lipid"):
            self.module.compare_dataset_sn(pd.concat([lipids, lipids.iloc[:1]]))

    def test_complete_lipids_excludes_missing_methods_and_invalid_sn(self):
        frame = pd.DataFrame([(dataset, lipid, 10)
                              for dataset in self.module.DATASETS for lipid in ["a", "b", "c"]],
                             columns=["Dataset", "Lipid", "S/N"])
        frame = frame.loc[~(frame.Dataset.eq("OTOT") & frame.Lipid.eq("b"))].copy()
        frame.loc[frame.Dataset.eq("OTIT") & frame.Lipid.eq("c"), "S/N"] = np.nan
        frame = pd.concat([frame, pd.DataFrame([(dataset, "d", 20) for dataset in self.module.DATASETS],
                                              columns=frame.columns)])
        complete = self.module.complete_lipids(frame)
        self.assertEqual(set(complete.Lipid), {"a", "d"})
        self.assertTrue(complete.groupby("Dataset").size().eq(2).all())

    def test_only_ns_brackets(self):
        from itertools import combinations
        lipids = self.module.prepare_lipid_sn(self.fixture())
        pairs = pd.DataFrame([{"Dataset 1": a, "Dataset 2": b,
                               "p": 1e-120, "Holm p": 1e-119, "Significant": True}
                              for a, b in combinations(self.module.DATASETS, 2)])
        pairs.loc[0, ["p", "Holm p", "Significant"]] = [0.03, 0.05, False]
        pairs.loc[9, ["p", "Holm p", "Significant"]] = [0.2, 0.4, False]
        anova = {"p": 1e-150, "Matched lipids": 180}
        fig = self.module.plot_dataset_sn(lipids, (anova, pairs))
        self.assertEqual(len(fig.axes[1].lines), 2)
        self.assertEqual(sum(text.get_text() == "ns" for text in fig.axes[1].texts), 2)

    def test_spectrum_first_intensity_and_sum(self):
        parse = self.module.spectrum_intensity
        self.assertEqual(parse("100.5:2 101.5:100 102.5:3", first_only=True), 2)
        self.assertEqual(parse("100.5:2 101.5:100 102.5:3"), 105)
        self.assertEqual(parse(" 100:1e3\n101:2.5e2 "), 1250)
        self.assertEqual(parse("100:0"), 0)
        self.assertEqual(parse("100:4 bad", first_only=True), 4)
        for invalid in [None, np.nan, "", "bad", "100:4 bad", "100:nan",
                        "100:inf", "100:-1", "nan:2", "100:2:3", "100:1e308 101:1e308"]:
            self.assertTrue(np.isnan(parse(invalid)), invalid)

    def test_intensity_filtering_and_lipid_medians(self):
        original = self.fixture()
        original["MS1 isotopic spectrum"] = "100:2 101:100"
        original["MS/MS spectrum"] = "100:2 101:100 102:3"
        original.loc[5, "MS1 isotopic spectrum"] = "100:10 101:999"
        original.loc[5, "MS/MS spectrum"] = "100:10 101:10"
        original.loc[6, "MS/MS spectrum"] = "100:0"
        for level, column, expected in [(1, "MS1 Intensity", 6), (2, "MS2 Summed Intensity", 62.5)]:
            lipids = self.module.prepare_lipid_intensity(original, level)
            value = lipids.loc[lipids.Dataset.eq("OTOT") & lipids.Lipid.eq("PC 34:1"), column].item()
            self.assertEqual(value, expected)
            self.assertEqual(lipids.loc[lipids.Lipid.eq("unique"), "Dataset"].tolist(), ["OTOT"])
            if level == 2:
                self.assertTrue(lipids.loc[lipids.Lipid.eq("unique"), column].isna().all())
            fig = self.module.plot_dataset_sn(lipids, value_column=column)
            self.assertEqual(fig.axes[0].get_ylabel(), f"{column} (log scale)")
        self.assertNotIn("MS1 Intensity", original)
        with self.assertRaisesRegex(ValueError, "MS1 isotopic spectrum"):
            self.module.prepare_lipid_intensity(self.fixture(), 1)


if __name__ == "__main__":
    unittest.main()
